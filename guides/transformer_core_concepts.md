<a name="top"></a>
# Building a Didactic Transformer from Scratch with PyTorch

This guide breaks down the process of building a Transformer model from scratch using PyTorch. We will construct each component piece by piece, with clear code and explanations that tie back to the foundational "Attention Is All You Need" paper.

This tutorial is designed for a learning environment like a Jupyter Notebook or Google Colab.

### Table of Contents
1.  [**Step 1: Input & Positional Embeddings**](#step-1-input-embeddings-and-positional-encoding)
    *   [A Deeper Look: The `d_model` Hyperparameter](#a-deeper-look-the-d_model-hyperparameter)
2.  [**Step 2: The Attention Mechanism**](#step-2-the-attention-mechanism)
3.  [**Step 3: Assembling the Encoder & Decoder Blocks**](#step-3-assembling-the-encoder--decoder-blocks)
4.  [**Step 4: Building the Full Transformer**](#step-4-building-the-full-transformer)

---

## Step 1: Input Embeddings and Positional Encoding
The first task in the Encoder is to convert input tokens into information-rich vectors. This involves two key components:

1.  **Input Embeddings:** Mapping each token (word or character) to a dense vector.
2.  **Positional Encoding:** Injecting information about the token's position in the sequence into that vector.

The goal is to create a final tensor that combines both semantic meaning and sequential order.

<details>
<summary><b>Click to view the PyTorch code for Input and Positional Embeddings</b></summary>

```python
import torch
import torch.nn as nn
import math

class InputEmbeddings(nn.Module):
    def __init__(self, d_model: int, vocab_size: int):
        """
        Initializes the Input Embedding layer.

        Args:
            d_model (int): The dimensionality of the model's embeddings (e.g., 512).
            vocab_size (int): The total number of unique tokens in the vocabulary.
        """
        super().__init__()
        self.d_model = d_model
        self.vocab_size = vocab_size
        self.embedding = nn.Embedding(vocab_size, d_model)

    def forward(self, x):
        """
        Forward pass for the embedding layer.
        The paper scales the embeddings by the square root of the model dimension.
        """
        return self.embedding(x) * math.sqrt(self.d_model)

class PositionalEncoding(nn.Module):
    def __init__(self, d_model: int, seq_len: int, dropout: float):
        """
        Initializes the Positional Encoding layer.

        Args:
            d_model (int): The dimensionality of the model.
            seq_len (int): The maximum sequence length the model can handle.
            dropout (float): The dropout rate.
        """
        super().__init__()
        self.d_model = d_model
        self.seq_len = seq_len
        self.dropout = nn.Dropout(dropout)

        # Create a matrix of shape (seq_len, d_model) to hold the positional encodings.
        pe = torch.zeros(seq_len, d_model)

        # Create a tensor representing the positions (0, 1, 2, ..., seq_len-1).
        position = torch.arange(0, seq_len, dtype=torch.float).unsqueeze(1) # (seq_len, 1)

        # Calculate the division term for the sine and cosine functions.
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * (-math.log(10000.0) / d_model))

        # Apply sine to even indices, cosine to odd indices.
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)

        # Add a batch dimension and register 'pe' as a buffer (not a trainable parameter).
        pe = pe.unsqueeze(0) # (1, seq_len, d_model)
        self.register_buffer('pe', pe)

    def forward(self, x):
        """
        Forward pass for Positional Encoding.
        Adds the positional encoding to the input embeddings up to the length of the input sequence.
        """
        x = x + (self.pe[:, :x.shape[1], :]).requires_grad_(False)
        return self.dropout(x)
```
</details>

### A Deeper Look: The `d_model` Hyperparameter
The `d_model` is a crucial hyperparameter that specifies the **size of the vector used to represent each token throughout the entire model.**

*   **Embedding Size:** It's the length of the vector the `nn.Embedding` layer outputs for each token.
*   **Hidden Size:** It's the dimensionality of the "working space" inside every subsequent layer.
*   **A Measure of Capacity:** You can think of `d_model` as the "bandwidth" or "richness" of information that can be stored for each token. A large `d_model` allows the model to encode more subtle and complex nuances. The original paper used `d_model = 512`.

> **Why is a large `d_model` essential for genomics?**
> In natural language, a word like "bank" has rich, pre-loaded meaning. In genomics, a token like `A` (Adenine) is almost meaningless by itself. Its **entire meaning is derived from its context**: its position in a codon, its proximity to a promoter, its role in a 3D fold, etc.
>
> Because the individual tokens have so little intrinsic meaning, the model needs an enormous `d_model` "workspace" for each token to accumulate all this necessary contextual information from across the sequence. Genomic foundation models like DNABERT and Geneformer use standard `d_model` values (e.g., 768 or 512) for this exact reason.

[Back to Top](#top)

---

## Step 2: The Attention Mechanism
The "magic" of the attention mechanism is that it allows every token to look at every other token and decide which ones are most important for understanding its own meaning. We will build this in three modular pieces.

1.  **`FeedForwardBlock`**: A simple helper network used within the main block.
2.  **`MultiHeadAttentionBlock`**: The main class that computes the attention scores.
3.  **`ResidualConnection`**: A helper that implements the "Add & Norm" steps, crucial for training deep networks.

<details>
<summary><b>Click to view the PyTorch code for Attention components</b></summary>

```python
import torch
import torch.nn as nn
import math

class FeedForwardBlock(nn.Module):
    def __init__(self, d_model: int, d_ff: int, dropout: float):
        """Initializes the Feed-Forward Block."""
        super().__init__()
        self.linear_1 = nn.Linear(d_model, d_ff) # W1 and b1
        self.dropout = nn.Dropout(dropout)
        self.linear_2 = nn.Linear(d_ff, d_model) # W2 and b2

    def forward(self, x):
        # (batch, seq_len, d_model) -> (batch, seq_len, d_ff) -> (batch, seq_len, d_model)
        return self.linear_2(self.dropout(torch.relu(self.linear_1(x))))

class MultiHeadAttentionBlock(nn.Module):
    def __init__(self, d_model: int, h: int, dropout: float):
        """Initializes the Multi-Head Attention Block."""
        super().__init__()
        self.d_model = d_model
        self.h = h
        assert d_model % h == 0, "d_model is not divisible by h"
        self.d_k = d_model // h
        self.w_q = nn.Linear(d_model, d_model)
        self.w_k = nn.Linear(d_model, d_model)
        self.w_v = nn.Linear(d_model, d_model)
        self.w_o = nn.Linear(d_model, d_model)
        self.dropout = nn.Dropout(dropout)

    @staticmethod
    def attention(query, key, value, mask, dropout: nn.Dropout):
        """The static attention calculation function."""
        d_k = query.shape[-1]
        attention_scores = (query @ key.transpose(-2, -1)) / math.sqrt(d_k)
        if mask is not None:
            attention_scores.masked_fill_(mask == 0, -1e9)
        attention_scores = attention_scores.softmax(dim=-1)
        if dropout is not None:
            attention_scores = dropout(attention_scores)
        return (attention_scores @ value), attention_scores

    def forward(self, q, k, v, mask):
        query = self.w_q(q)
        key = self.w_k(k)
        value = self.w_v(v)
        query = query.view(query.shape[0], query.shape[1], self.h, self.d_k).transpose(1, 2)
        key = key.view(key.shape[0], key.shape[1], self.h, self.d_k).transpose(1, 2)
        value = value.view(value.shape[0], value.shape[1], self.h, self.d_k).transpose(1, 2)
        x, self.attention_scores = MultiHeadAttentionBlock.attention(query, key, value, mask, self.dropout)
        x = x.transpose(1, 2).contiguous().view(x.shape[0], -1, self.h * self.d_k)
        return self.w_o(x)

class LayerNormalization(nn.Module):
    def __init__(self, eps: float = 10**-6):
        """Initializes a custom Layer Normalization module for didactic clarity."""
        super().__init__()
        self.eps = eps
        self.alpha = nn.Parameter(torch.ones(1)) # Learnable scale
        self.bias = nn.Parameter(torch.zeros(1))  # Learnable shift

    def forward(self, x):
        mean = x.mean(dim=-1, keepdim=True)
        std = x.std(dim=-1, keepdim=True)
        return self.alpha * (x - mean) / (std + self.eps) + self.bias

class ResidualConnection(nn.Module):
    def __init__(self, dropout: float):
        """Initializes the Residual Connection (Add & Norm) module."""
        super().__init__()
        self.dropout = nn.Dropout(dropout)
        self.norm = LayerNormalization()

    def forward(self, x, sublayer):
        """Applies the residual connection to a sublayer.
        This is a Pre-LN style: x + Dropout(Sublayer(Norm(x)))
        """
        return x + self.dropout(sublayer(self.norm(x)))
```
</details>

[Back to Top](#top)

---

## Step 3: Assembling the Encoder & Decoder Blocks
Now we compose the pieces from Step 2 into complete Encoder and Decoder blocks.

An **Encoder Block** has two sub-layers:
1.  A Multi-Head Attention block.
2.  A Feed-Forward block.

A **Decoder Block** is similar but has three sub-layers:
1.  A "Masked" Multi-Head Self-Attention block (looks at what it has generated so far).
2.  A Multi-Head Cross-Attention block (looks at the Encoder's output).
3.  A Feed-Forward block.

<details>
<summary><b>Click to view the PyTorch code for Encoder and Decoder Blocks</b></summary>

```python
import torch.nn as nn

# Assumes MultiHeadAttentionBlock, FeedForwardBlock, ResidualConnection are defined

class EncoderBlock(nn.Module):
    def __init__(self, self_attention_block: MultiHeadAttentionBlock, feed_forward_block: FeedForwardBlock, dropout: float):
        """Initializes a single Encoder Block."""
        super().__init__()
        self.self_attention_block = self_attention_block
        self.feed_forward_block = feed_forward_block
        self.residual_connections = nn.ModuleList([ResidualConnection(dropout) for _ in range(2)])

    def forward(self, x, src_mask):
        # In self-attention, query, key, and value are all the same 'x'.
        x = self.residual_connections[0](x, lambda x: self.self_attention_block(x, x, x, src_mask))
        x = self.residual_connections[1](x, self.feed_forward_block)
        return x

class DecoderBlock(nn.Module):
    def __init__(self, self_attention_block: MultiHeadAttentionBlock, cross_attention_block: MultiHeadAttentionBlock, feed_forward_block: FeedForwardBlock, dropout: float):
        """Initializes a single Decoder Block."""
        super().__init__()
        self.self_attention_block = self_attention_block
        self.cross_attention_block = cross_attention_block
        self.feed_forward_block = feed_forward_block
        self.residual_connections = nn.ModuleList([ResidualConnection(dropout) for _ in range(3)])

    def forward(self, x, encoder_output, src_mask, tgt_mask):
        # 1. Masked Self-Attention (looks at itself)
        x = self.residual_connections[0](x, lambda x: self.self_attention_block(x, x, x, tgt_mask))
        # 2. Cross-Attention (Query from decoder, Key/Value from encoder)
        x = self.residual_connections[1](x, lambda x: self.cross_attention_block(x, encoder_output, encoder_output, src_mask))
        # 3. Feed-Forward Network
        x = self.residual_connections[2](x, self.feed_forward_block)
        return x
```
</details>

[Back to Top](#top)

---

## Step 4: Building the Full Transformer
Finally, we assemble everything into the main `Transformer` class. This involves:
1.  **`Encoder`**: A class to stack multiple `EncoderBlock` layers.
2.  **`Decoder`**: A class to stack multiple `DecoderBlock` layers.
3.  **`ProjectionLayer`**: The final layer to map the output vectors to word probabilities.
4.  **`Transformer`**: The main class that orchestrates the entire process.

<details>
<summary><b>Click to view the PyTorch code for the full Transformer assembly</b></summary>

```python
import torch
import torch.nn as nn

# Assumes all previous component classes are defined

class Encoder(nn.Module):
    def __init__(self, layers: nn.ModuleList):
        """Initializes the full Encoder stack."""
        super().__init__()
        self.layers = layers
        self.norm = LayerNormalization()

    def forward(self, x, mask):
        for layer in self.layers:
            x = layer(x, mask)
        return self.norm(x)

class Decoder(nn.Module):
    def __init__(self, layers: nn.ModuleList):
        """Initializes the full Decoder stack."""
        super().__init__()
        self.layers = layers
        self.norm = LayerNormalization()

    def forward(self, x, encoder_output, src_mask, tgt_mask):
        for layer in self.layers:
            x = layer(x, encoder_output, src_mask, tgt_mask)
        return self.norm(x)

class ProjectionLayer(nn.Module):
    def __init__(self, d_model: int, vocab_size: int):
        """Initializes the final projection layer to vocabulary size."""
        super().__init__()
        self.proj = nn.Linear(d_model, vocab_size)

    def forward(self, x):
        # log_softmax is often used for better numerical stability during training.
        return torch.log_softmax(self.proj(x), dim=-1)

class Transformer(nn.Module):
    def __init__(self, encoder: Encoder, decoder: Decoder, src_embed: InputEmbeddings, tgt_embed: InputEmbeddings, src_pos: PositionalEncoding, tgt_pos: PositionalEncoding, projection_layer: ProjectionLayer):
        """Initializes the full Transformer model."""
        super().__init__()
        self.encoder = encoder
        self.decoder = decoder
        self.src_embed = src_embed
        self.tgt_embed = tgt_embed
        self.src_pos = src_pos
        self.tgt_pos = tgt_pos
        self.projection_layer = projection_layer

    def encode(self, src, src_mask):
        """Performs the encoding part of the model."""
        src = self.src_embed(src)
        src = self.src_pos(src)
        return self.encoder(src, src_mask)

    def decode(self, encoder_output, src_mask, tgt, tgt_mask):
        """Performs the decoding part of the model."""
        tgt = self.tgt_embed(tgt)
        tgt = self.tgt_pos(tgt)
        return self.decoder(tgt, encoder_output, src_mask, tgt_mask)

    def project(self, x):
        """Applies the final projection layer."""
        return self.projection_layer(x)
```
</details>

With these final classes, we have successfully built a complete Transformer from the ground up. Every piece from the original paper's diagram is accounted for, from the initial embedding to the final prediction.

[Back to Top](#top)

## Credits
This lesson was provided by Gemini Pro (AI).
