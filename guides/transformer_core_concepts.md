### **Project: Building a Didactic Transformer from Scratch**

Our goal is to write a simple, clean, and well-commented Transformer in Python using the PyTorch library. We will build it piece by piece, and at each step, I will provide the code and explain exactly how it corresponds to the theory we've discussed.

**The Environment:**
For this, a simple environment like a Jupyter Notebook or Google Colab is perfect. Google Colab is often easiest as it requires no local setup and has PyTorch pre-installed.

**Let's begin with the very first step in the Encoder:** taking the input tokens and converting them into information-rich vectors. This involves two key components we've discussed:

1.  **Input Embeddings:** Mapping each token (word/sub-word) to a dense vector.
2.  **Positional Encoding:** Adding information about the token's position in the sequence to that vector.

Here is the code for a PyTorch module that accomplishes both. I've included detailed comments to explain each line.

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
        # The main embedding layer that maps token IDs to vectors.
        self.embedding = nn.Embedding(vocab_size, d_model)

    def forward(self, x):
        """
        Forward pass for the embedding layer.

        Args:
            x (torch.Tensor): A tensor of token IDs with shape (batch_size, seq_len).

        Returns:
            torch.Tensor: The embedded vectors, scaled by sqrt(d_model).
                          Shape: (batch_size, seq_len, d_model).
        """
        # The paper scales the embeddings by the square root of the model dimension.
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
        # This is the `10000^(2i/d_model)` term from the paper.
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * (-math.log(10000.0) / d_model))

        # Apply the sine function to even indices of the 'pe' matrix.
        pe[:, 0::2] = torch.sin(position * div_term)
        # Apply the cosine function to odd indices of the 'pe' matrix.
        pe[:, 1::2] = torch.cos(position * div_term)

        # Add a batch dimension to the positional encoding matrix so it can be added to the input embeddings.
        pe = pe.unsqueeze(0) # (1, seq_len, d_model)

        # Register 'pe' as a buffer. A buffer is part of the model's state, but it is not a parameter
        # to be trained. This is important because positional encodings are fixed.
        self.register_buffer('pe', pe)

    def forward(self, x):
        """
        Forward pass for Positional Encoding.

        Args:
            x (torch.Tensor): The input embeddings with shape (batch_size, seq_len, d_model).

        Returns:
            torch.Tensor: The embeddings with positional information added.
                          Shape: (batch_size, seq_len, d_model).
        """
        # Add the positional encoding to the input embeddings.
        # x.shape[1] is the sequence length of the current batch, which might be smaller than self.seq_len.
        # We add the positional encodings up to the length of the input sequence.
        x = x + (self.pe[:, :x.shape[1], :]).requires_grad_(False)
        return self.dropout(x)

```

### The Goal: Addition Requires Identical Shapes

The entire purpose of the `PositionalEncoding` module is to create a tensor that has the **exact same shape** as the output of the `InputEmbeddings` module, so that they can be added together element-wise.

Let's define our key dimensions first. These are **hyperparameters**—design choices you make before training begins. They are the "prior values" you mentioned.

*   `batch_size`: The number of sentences we process at once (e.g., 32).
*   `seq_len`: The length of the sequences in the current batch (e.g., 100 tokens).
*   `d_model`: The embedding dimension of our model (e.g., 512). This is a core hyperparameter.
*   `max_seq_len`: The absolute maximum sequence length our model is designed to handle (e.g., 5000).


### Dimensionality Dynamics: A Step-by-Step Walkthrough

Let's trace a tensor `x` through the process.

**Step 1: Raw Input**
The input `x` is a batch of sentences where each word has been converted to its token ID.
*   **Shape of `x`:** `(batch_size, seq_len)`
*   **Example:** `(32, 100)` -> A batch of 32 sentences, each 100 tokens long.

**Step 2: Inside `InputEmbeddings.forward(x)`**
The `nn.Embedding` layer is essentially a big lookup table. For each of the `100` token IDs in a sentence, it looks up its corresponding vector of size `d_model`.
*   **Operation:** `self.embedding(x)`
*   **Shape after `embedding`:** `(batch_size, seq_len, d_model)`
*   **Example:** `(32, 100, 512)` -> We now have 32 sentences, each composed of 100 vectors of size 512.
*   The multiplication by `math.sqrt(self.d_model)` is a scalar operation, so it doesn't change the shape.

**This shape, `(batch_size, seq_len, d_model)`, is our target shape.** The positional encoding must match this.

---

**Step 3: Inside `PositionalEncoding.__init__` (The One-Time Setup)**
This part is tricky because it's creating a large, reusable matrix that will work for *any* batch.

*   **`pe = torch.zeros(max_seq_len, d_model)`**
    *   **Why?** We create a big empty matrix to hold encodings for the longest possible sentence we might ever see.
    *   **Shape:** `(5000, 512)`

*   **`position = torch.arange(...).unsqueeze(1)`**
    *   `arange` creates `[0, 1, 2, ..., 4999]`. Shape: `(5000)`.
    *   `unsqueeze(1)` adds a dimension. Shape: `(5000, 1)`.
    *   **Why?** This is crucial for **broadcasting**. It allows us to multiply this `(5000, 1)` position vector with the `(256)` frequency vector to produce a final `(5000, 256)` matrix without a loop.

*   **`div_term = ...`**
    *   This creates the frequencies for the sine/cosine waves.
    *   **Shape:** `(d_model / 2)` -> `(256)`

*   **`pe[:, 0::2] = torch.sin(position * div_term)`**
    *   This is where the magic happens. `position` `(5000, 1)` is multiplied with `div_term` `(256)`. Thanks to broadcasting, the result is a matrix of shape `(5000, 256)`. This fills all the even columns of `pe`.
    *   The `cos` part does the same for the odd columns.

*   **`pe = pe.unsqueeze(0)`**
    *   This is the final critical setup step. We add a batch dimension to our `pe` matrix.
    *   **Why?** So its dimensionality aligns with our input tensor `x`, which also has a batch dimension.
    *   **Final Shape of `self.pe` buffer:** `(1, max_seq_len, d_model)` -> `(1, 5000, 512)`

---

**Step 4: Inside `PositionalEncoding.forward(x)` (The Final Addition)**
The embedded input `x` arrives from Step 2.

*   **Shape of input `x`:** `(batch_size, seq_len, d_model)` -> `(32, 100, 512)`

Now, we perform the operation `x + self.pe[:, :x.shape[1], :]`. Let's break down the **slicing** on the right side:

*   `self.pe` has shape `(1, 5000, 512)`.
*   `self.pe[:, :x.shape[1], :]` means:
    *   `:` in the 1st dim: Keep the batch dimension.
    *   `:x.shape[1]` in the 2nd dim: Slice the sequence from `0` to `100` (the length of our current input).
    *   `:` in the 3rd dim: Keep the full `d_model` dimension.
*   **Shape of the slice:** `(1, 100, 512)`

Now we can do the addition:
`x` `(32, 100, 512)` + `pe_slice` `(1, 100, 512)`

**Why this works:** This is the final piece of broadcasting magic. PyTorch sees the `1` in the batch dimension of `pe_slice` and automatically "stretches" or "clones" it 32 times to match the batch dimension of `x`.

*   **Shape of the final output:** `(32, 100, 512)`

We have successfully added the positional information without changing the shape of our data tensor.

### Summary Table

| Step | Operation | Output Tensor Shape | Why? |
| :--- | :--- | :--- | :--- |
| 1 | Raw Input `x` | `(batch, seq_len)` | Batch of token IDs |
| 2 | `nn.Embedding` | `(batch, seq_len, d_model)` | Convert IDs to dense vectors |
| 3 | `PositionalEncoding` setup | `(1, max_seq_len, d_model)` | Create a reusable buffer for all positions |
| 4 | **The `forward` pass** | | |
| | Input `x` from Step 2 | `(batch, seq_len, d_model)` | |
| | Slice `pe` buffer | `(1, seq_len, d_model)` | Match the current batch's sequence length |
| | Add `x + pe_slice` | `(batch, seq_len, d_model)` | Broadcasting adds the `(1, ...)` slice to every item in the batch |

This detailed tracking shows that the dimensions are not arbitrary. They are a careful chain of transformations designed to make the final, crucial addition operation possible, using **hyperparameters** as the fixed anchors and **broadcasting/slicing** as the dynamic tools.

`d_model` **(dimensionality of the model)** is a hyperparameter that specifies the **size of the vector used to represent each token throughout the entire model.**

*   **Embedding Size:** It's the length of the vector that the `nn.Embedding` layer outputs for each token.
*   **Hidden Size:** It's the dimensionality of the "working space" inside every subsequent layer of the Transformer. The `Query`, `Key`, and `Value` vectors will all have a size related to `d_model`. The output of each attention layer and each feed-forward network will be a vector of size `d_model`.
*   **A Measure of "Capacity":** You can think of `d_model` as the "bandwidth" or "richness" of the information that can be stored for each token.
    *   A small `d_model` (e.g., 64) can only capture coarse relationships.
    *   A large `d_model` (e.g., 512, 1024, or even larger in huge models) allows the model to encode much more subtle and complex nuances about a token's meaning, its grammatical role, and its relationship to other tokens.

**The Trade-off:**
*   **Larger `d_model`:** Leads to a more powerful and expressive model, but...
    *   It dramatically increases the number of parameters (model size).
    *   It increases the computational cost (memory and processing time).
    *   It requires more data to train effectively without overfitting.

In the original "Attention Is All You Need" paper, the value used was **`d_model = 512`**. This has become a common baseline for "base" sized models. "Large" models often use `d_model = 1024` or higher.

### The Vocabulary vs. The Context

In natural language, a word like "bank" has a rich, pre-loaded meaning. The `nn.Embedding` layer learns a unique vector for "bank" that is already quite informative. The rest of the model's job is to figure out whether it's a river bank or a financial bank.

In genomics, you have a tiny vocabulary: `A`, `C`, `G`, `T`, (and perhaps `N` for unknown).
*   A token `A` (Adenine) by itself has very little meaning. It's just a letter.
*   The **entire meaning** of that `A` is derived from its context: the vast, complex sequence of millions of other bases around it. Its meaning is defined by its position within a codon, its proximity to a promoter region, its role in a gene, its 3D folding implications, etc.

### Why a Large `d_model` is Crucial for DNA

Because the individual tokens are almost meaningless, the model needs an enormous "workspace" or "bandwidth" (`d_model`) for each token to accumulate all of this necessary contextual information.

1.  **Capturing Long-Range Dependencies:** A single `A` at position 1,000 might be functionally linked to a `G` at position 500,000 due to how the DNA folds in 3D space. The model needs a large `d_model` vector at each position to hold the complex "summary" of these incredibly long-range interactions. A small vector would quickly become saturated and unable to store more information.

2.  **Encoding Positional and Structural Information:** The `d_model` vector for a given `T` doesn't just need to encode "this is Thymine." It needs to encode:
    *   "This is a `T` at the third position of a codon."
    *   "This `T` is part of a TATA-box, a key promoter region."
    *   "This `T` is in a region that is highly conserved across species."
    *   "This `T`, if mutated to a `C`, is associated with a specific disease."
    A large `d_model` provides the capacity to store all these different layers of annotation simultaneously.

3.  **Compensating for a "Dumb" Embedding Layer:** The `nn.Embedding` layer for DNA is very small. It only has to learn 4-5 unique vectors. It's not doing much work. All the "heavy lifting" of creating meaning must be done by the subsequent Transformer layers (attention and feed-forward networks), and these layers operate on the `d_model` dimension.

### The Evidence in Practice

Specialized genomic foundation models like **DNABERT**, **Geneformer**, and others use the exact same Transformer architecture. And yes, they use standard `d_model` values:
*   DNABERT uses `d_model = 768` (the same as BERT-base).
*   Geneformer uses `d_model = 512`.

They have tiny vocabularies, but they are processing sequences that can be tens or hundreds of thousands of tokens long, where the context is everything. The large `d_model` is not for representing the token itself, but for representing the **token-in-its-unfathomably-complex-context.**

The "magic" of the attention mechanism is that it allows every token to look at every other token in the sequence and decide which ones are most important for understanding its own meaning. The "Multi-Head" part is a clever trick to do this from multiple different perspectives at the same time.

We will build this in three logical, modular pieces, just as a software engineer would:
1.  **`FeedForwardBlock`**: A simple helper network that is used within the attention block.
2.  **`MultiHeadAttentionBlock`**: The main class that computes the attention scores.
3.  **`ResidualConnection`**: Another helper that implements the "add & norm" steps from the paper's diagram.

---

### Code Part 1: The Feed-Forward Block

After the model calculates attention, it passes the result through a simple two-layer neural network. This is done for each token's vector individually. This step adds more learnable parameters and allows the model to process the attention output in a more complex way.

```python
import torch.nn as nn

class FeedForwardBlock(nn.Module):
    def __init__(self, d_model: int, d_ff: int, dropout: float):
        """
        Initializes the Feed-Forward Block.

        Args:
            d_model (int): The dimensionality of the model.
            d_ff (int): The dimensionality of the inner feed-forward layer. A common value is 4 * d_model.
            dropout (float): The dropout rate.
        """
        super().__init__()
        self.linear_1 = nn.Linear(d_model, d_ff) # W1 and b1
        self.dropout = nn.Dropout(dropout)
        self.linear_2 = nn.Linear(d_ff, d_model) # W2 and b2

    def forward(self, x):
        """
        Forward pass for the Feed-Forward Block.

        Input shape: (batch_size, seq_len, d_model)
        Output shape: (batch_size, seq_len, d_model)
        """
        # (batch, seq_len, d_model) --> (batch, seq_len, d_ff) --> (batch, seq_len, d_model)
        return self.linear_2(self.dropout(torch.relu(self.linear_1(x))))

```

### Code Part 2: The Multi-Head Attention Block

This is the main event. It takes the input vectors and creates the `Query`, `Key`, and `Value` matrices. It then calculates the attention scores and the final output.

```python
import torch
import torch.nn as nn
import math

class MultiHeadAttentionBlock(nn.Module):
    def __init__(self, d_model: int, h: int, dropout: float):
        """
        Initializes the Multi-Head Attention Block.

        Args:
            d_model (int): The dimensionality of the model.
            h (int): The number of attention heads.
            dropout (float): The dropout rate.
        """
        super().__init__()
        self.d_model = d_model
        self.h = h
        # Ensure that the model dimension is divisible by the number of heads
        assert d_model % h == 0, "d_model is not divisible by h"

        # d_k is the dimension of each head's key/query/value vectors
        self.d_k = d_model // h

        # The linear layers for Query, Key, Value, and the final output
        self.w_q = nn.Linear(d_model, d_model) # Wq
        self.w_k = nn.Linear(d_model, d_model) # Wk
        self.w_v = nn.Linear(d_model, d_model) # Wv
        self.w_o = nn.Linear(d_model, d_model) # Wo

        self.dropout = nn.Dropout(dropout)

    @staticmethod
    def attention(query, key, value, mask, dropout: nn.Dropout):
        """
        The static attention calculation function.

        Args:
            query, key, value: Tensors of shape (batch, h, seq_len, d_k)
            mask: A mask to hide certain interactions (e.g., padding or future tokens).
            dropout: The dropout layer.

        Returns:
            A tuple of (output tensor, attention scores tensor).
        """
        d_k = query.shape[-1]

        # (batch, h, seq_len, d_k) @ (batch, h, d_k, seq_len) --> (batch, h, seq_len, seq_len)
        attention_scores = (query @ key.transpose(-2, -1)) / math.sqrt(d_k)

        # Apply mask before the softmax
        if mask is not None:
            attention_scores.masked_fill_(mask == 0, -1e9)

        attention_scores = attention_scores.softmax(dim=-1) # (batch, h, seq_len, seq_len)

        if dropout is not None:
            attention_scores = dropout(attention_scores)

        # (batch, h, seq_len, seq_len) @ (batch, h, seq_len, d_k) --> (batch, h, seq_len, d_k)
        return (attention_scores @ value), attention_scores

    def forward(self, q, k, v, mask):
        """
        Forward pass for the Multi-Head Attention Block.

        Input shapes: (batch_size, seq_len, d_model) for q, k, v
        Output shape: (batch_size, seq_len, d_model)
        """
        query = self.w_q(q) # (batch, seq_len, d_model)
        key = self.w_k(k)   # (batch, seq_len, d_model)
        value = self.w_v(v) # (batch, seq_len, d_model)

        # Reshape for multi-head processing:
        # (batch, seq_len, d_model) --> (batch, seq_len, h, d_k) --> (batch, h, seq_len, d_k)
        query = query.view(query.shape[0], query.shape[1], self.h, self.d_k).transpose(1, 2)
        key = key.view(key.shape[0], key.shape[1], self.h, self.d_k).transpose(1, 2)
        value = value.view(value.shape[0], value.shape[1], self.h, self.d_k).transpose(1, 2)

        # Calculate attention
        x, self.attention_scores = MultiHeadAttentionBlock.attention(query, key, value, mask, self.dropout)

        # Combine the heads back together
        # (batch, h, seq_len, d_k) --> (batch, seq_len, h, d_k) --> (batch, seq_len, d_model)
        x = x.transpose(1, 2).contiguous().view(x.shape[0], -1, self.h * self.d_k)

        # Final linear layer
        # (batch, seq_len, d_model) --> (batch, seq_len, d_model)
        return self.w_o(x)
```

### Code Part 3: The Residual Connection

This component implements the "Add & Norm" step shown in the diagrams. It takes the output of a sub-layer (like Multi-Head Attention), adds the original input back to it (a "residual" or "skip" connection), and then normalizes the result. This is crucial for training deep networks.

```python
import torch.nn as nn

class LayerNormalization(nn.Module):
    def __init__(self, eps: float = 10**-6):
        """
        Initializes the Layer Normalization module.

        Args:
            eps (float): A small value added for numerical stability.
        """
        super().__init__()
        self.eps = eps
        self.alpha = nn.Parameter(torch.ones(1)) # Learnable scale parameter
        self.bias = nn.Parameter(torch.zeros(1))  # Learnable shift parameter

    def forward(self, x):
        mean = x.mean(dim = -1, keepdim=True)
        std = x.std(dim = -1, keepdim=True)
        return self.alpha * (x - mean) / (std + self.eps) + self.bias

class ResidualConnection(nn.Module):
    def __init__(self, dropout: float):
        """
        Initializes the Residual Connection module with Layer Normalization.
        """
        super().__init__()
        self.dropout = nn.Dropout(dropout)
        self.norm = LayerNormalization()

    def forward(self, x, sublayer):
        """
        Forward pass for the residual connection.

        Args:
            x: The original input to the sublayer.
            sublayer: The sublayer function (e.g., the multi-head attention block).

        Returns:
            The output after the residual connection and layer normalization.
        """
        # The architecture is x + Dropout(Sublayer(Norm(x)))
        # This is a slight variation from the original paper, often called "Pre-LN"
        # and found to be more stable for training deep transformers.
        return x + self.dropout(sublayer(self.norm(x)))
```
*(Note: I've included a custom `LayerNormalization` class for didactic clarity, though in practice one would often use `nn.LayerNorm`)*.

Do you see how the `MultiHeadAttentionBlock.forward` function orchestrates the creation of Q, K, and V? And how the `attention` static method performs the core `(Q * K^T) / sqrt(d_k)` calculation we discussed in theory?


### A Direct Comparison: Our PyTorch vs. Keras

Here is how you would implement the exact same Multi-Head Attention block we just built, but using the high-level Keras API.

**Our PyTorch `MultiHeadAttentionBlock` (Low-Level):**

```python
# Our detailed, manual implementation
class MultiHeadAttentionBlock(nn.Module):
    def __init__(self, d_model: int, h: int, dropout: float):
        super().__init__()
        # ... all the manual setup ...
        self.d_k = d_model // h
        self.w_q = nn.Linear(d_model, d_model)
        self.w_k = nn.Linear(d_model, d_model)
        self.w_v = nn.Linear(d_model, d_model)
        self.w_o = nn.Linear(d_model, d_model)

    def forward(self, q, k, v, mask):
        # ... manually project q, k, v ...
        query = self.w_q(q)
        # ... manually reshape for multi-head ...
        query = query.view(query.shape[0], query.shape[1], self.h, self.d_k).transpose(1, 2)
        # ... (repeat for key and value) ...
        # ... manually call the attention function ...
        x, self.attention_scores = MultiHeadAttentionBlock.attention(query, key, value, mask, self.dropout)
        # ... manually reshape back ...
        x = x.transpose(1, 2).contiguous().view(x.shape[0], -1, self.h * self.d_k)
        # ... manually apply final linear layer ...
        return self.w_o(x)

```

**The Keras `MultiHeadAttention` Layer (High-Level):**

```python
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers

# Define hyperparameters
d_model = 512
num_heads = 8

# Create the layer
# The layer itself handles creating W_q, W_k, W_v, and W_o internally.
# It knows how to split d_model into num_heads * key_dim.
mha_layer = layers.MultiHeadAttention(
    num_heads=num_heads,
    key_dim=d_model // num_heads,  # This is our d_k
    output_shape=d_model
)

# Assume we have input tensors 'query', 'key', 'value'
# Shape: (batch_size, seq_len, d_model)
# You just call the layer like a function.
attention_output = mha_layer(
    query=query_input,
    key=key_input,
    value=value_input,
    attention_mask=mask # The mask is passed directly
)

# 'attention_output' now has the correct shape: (batch_size, seq_len, d_model)
```

### The Keras Perspective: Key Abstractions

1.  **No Manual Reshaping:** Notice the complete absence of `.view()`, `.transpose()`, or `.contiguous()`. You provide tensors of shape `(batch, seq, d_model)`, and Keras handles all the internal splitting into heads, calculating attention, and concatenating the results.

2.  **Internal Weight Creation:** You don't define `w_q`, `w_k`, `w_v`, or `w_o`. The `MultiHeadAttention` layer knows it needs these and creates them automatically when you initialize it. It's a self-contained "black box" that does one job.

3.  **Focus on Intent, Not Mechanics:** The Keras code reads like a description of what you *want* to do ("I want to perform multi-head attention"), not a step-by-step instruction manual of *how* to do it.

### Code Part 4: Assembling the Encoder Layer

The diagram for an Encoder Layer shows two main sub-layers:
1.  A Multi-Head Attention block.
2.  A Feed-Forward block.
Each of these is wrapped in a "Residual Connection" (the "Add & Norm" step).

Let's write the code for one complete Encoder block. Notice how it's simply composing the pieces we've already built.

```python
import torch.nn as nn

class EncoderBlock(nn.Module):
    def __init__(self, self_attention_block: MultiHeadAttentionBlock, feed_forward_block: FeedForwardBlock, dropout: float):
        """
        Initializes a single Encoder Block.

        Args:
            self_attention_block (MultiHeadAttentionBlock): The multi-head self-attention mechanism.
            feed_forward_block (FeedForwardBlock): The feed-forward network.
            dropout (float): The dropout rate.
        """
        super().__init__()
        self.self_attention_block = self_attention_block
        self.feed_forward_block = feed_forward_block
        # We need two residual connections, one for each sub-layer.
        self.residual_connections = nn.ModuleList([ResidualConnection(dropout) for _ in range(2)])

    def forward(self, x, src_mask):
        """
        Forward pass for the Encoder Block.

        Args:
            x (torch.Tensor): The input tensor from the previous layer.
            src_mask: The mask to hide padding tokens in the source sequence.

        Returns:
            The output tensor of the same shape as the input.
        """
        # First residual connection: applies self-attention
        # Note: In self-attention, query, key, and value are all the same 'x'.
        x = self.residual_connections[0](x, lambda x: self.self_attention_block(x, x, x, src_mask))

        # Second residual connection: applies the feed-forward network
        x = self.residual_connections[1](x, self.feed_forward_block)

        return x
```
**Clarity Note:** The `lambda x: ...` syntax is just a clean way to pass our `self_attention_block` as a function to the `ResidualConnection` module, ensuring `x` is used for `q`, `k`, and `v`.

---

### Code Part 5: Assembling the Decoder Layer

The Decoder Layer is very similar, but with one key difference: it has **three** sub-layers.
1.  A "Masked" Multi-Head Self-Attention block (to look at what it has written so far).
2.  A Multi-Head Cross-Attention block (to look at the Encoder's output).
3.  A Feed-Forward block.

Each of these three is also wrapped in a residual connection.

```python
import torch.nn as nn

class DecoderBlock(nn.Module):
    def __init__(self, self_attention_block: MultiHeadAttentionBlock, cross_attention_block: MultiHeadAttentionBlock, feed_forward_block: FeedForwardBlock, dropout: float):
        """
        Initializes a single Decoder Block.

        Args:
            self_attention_block (MultiHeadAttentionBlock): The self-attention for the target sequence.
            cross_attention_block (MultiHeadAttentionBlock): The cross-attention to look at the encoder output.
            feed_forward_block (FeedForwardBlock): The feed-forward network.
            dropout (float): The dropout rate.
        """
        super().__init__()
        self.self_attention_block = self_attention_block
        self.cross_attention_block = cross_attention_block
        self.feed_forward_block = feed_forward_block
        # We now need three residual connections.
        self.residual_connections = nn.ModuleList([ResidualConnection(dropout) for _ in range(3)])

    def forward(self, x, encoder_output, src_mask, tgt_mask):
        """
        Forward pass for the Decoder Block.

        Args:
            x (torch.Tensor): The input from the target sequence (from previous decoder layer).
            encoder_output (torch.Tensor): The output from the top of the encoder stack.
            src_mask: The mask for the encoder output.
            tgt_mask: The mask for the target sequence (hides padding AND future tokens).

        Returns:
            The output tensor of the same shape as the input 'x'.
        """
        # 1. Masked Self-Attention (looks at itself)
        # Query, Key, Value are all 'x' from the target sequence. tgt_mask is used.
        x = self.residual_connections[0](x, lambda x: self.self_attention_block(x, x, x, tgt_mask))

        # 2. Cross-Attention (looks at the encoder's output)
        # Query comes from the decoder ('x'), but Key and Value come from the encoder_output.
        # The src_mask is used here to hide padding in the encoder's output.
        x = self.residual_connections[1](x, lambda x: self.cross_attention_block(x, encoder_output, encoder_output, src_mask))

        # 3. Feed-Forward Network
        x = self.residual_connections[2](x, self.feed_forward_block)

        return x
```
We will create several final classes:

1.  **`Encoder`**: A class that simply stacks multiple `EncoderBlock` layers on top of each other.
2.  **`Decoder`**: A class that stacks multiple `DecoderBlock` layers.
3.  **`ProjectionLayer`**: The final layer of the model that converts the decoder's output vectors back into word probabilities.
4.  **`Transformer`**: The main class that orchestrates the entire process, containing an Encoder, a Decoder, and all the necessary embedding and projection layers.

---

### Code Part 6: The Full Encoder and Decoder Stacks

These classes are quite simple. Their main job is to hold a list of our `EncoderBlock` or `DecoderBlock` modules and pass the data through them sequentially.

```python
import torch.nn as nn

class Encoder(nn.Module):
    def __init__(self, layers: nn.ModuleList):
        """
        Initializes the full Encoder stack.

        Args:
            layers (nn.ModuleList): A list of EncoderBlock layers.
        """
        super().__init__()
        self.layers = layers
        # The final normalization layer for the encoder stack
        self.norm = LayerNormalization()

    def forward(self, x, mask):
        # Pass the input through each layer in the stack
        for layer in self.layers:
            x = layer(x, mask)
        # Apply final normalization
        return self.norm(x)

class Decoder(nn.Module):
    def __init__(self, layers: nn.ModuleList):
        """
        Initializes the full Decoder stack.

        Args:
            layers (nn.ModuleList): A list of DecoderBlock layers.
        """
        super().__init__()
        self.layers = layers
        # The final normalization layer for the decoder stack
        self.norm = LayerNormalization()

    def forward(self, x, encoder_output, src_mask, tgt_mask):
        # Pass the input through each layer in the stack
        for layer in self.layers:
            x = layer(x, encoder_output, src_mask, tgt_mask)
        # Apply final normalization
        return self.norm(x)
```

### Code Part 7: The Final Projection Layer

After the decoder stack produces its final output vectors (shape: `(batch, seq_len, d_model)`), we need to convert these back into vocabulary words. This is done by a final linear layer followed by a softmax function.

```python
import torch.nn as nn

class ProjectionLayer(nn.Module):
    def __init__(self, d_model: int, vocab_size: int):
        """
        Initializes the final projection layer.

        Args:
            d_model (int): The dimensionality of the model.
            vocab_size (int): The size of the target vocabulary.
        """
        super().__init__()
        # A linear layer to project from d_model to the vocabulary size
        self.proj = nn.Linear(d_model, vocab_size)

    def forward(self, x):
        """
        Forward pass for the projection layer.

        Input shape: (batch_size, seq_len, d_model)
        Output shape: (batch_size, seq_len, vocab_size) - These are called logits.
        """
        # The log_softmax is often used for better numerical stability during training.
        return torch.log_softmax(self.proj(x), dim=-1)
```

---

### Code Part 8: The Grand Finale - The Full Transformer

Now we assemble everything into our main `Transformer` class. This class will contain all the components we have built, from the initial embeddings to the final projection layer, and will define the end-to-end forward pass.

```python
import torch
import torch.nn as nn

class Transformer(nn.Module):
    def __init__(self, encoder: Encoder, decoder: Decoder, src_embed: InputEmbeddings, tgt_embed: InputEmbeddings, src_pos: PositionalEncoding, tgt_pos: PositionalEncoding, projection_layer: ProjectionLayer):
        """
        Initializes the full Transformer model.
        """
        super().__init__()
        self.encoder = encoder
        self.decoder = decoder
        self.src_embed = src_embed # Embeddings for the source language
        self.tgt_embed = tgt_embed # Embeddings for the target language
        self.src_pos = src_pos     # Positional encoding for the source
        self.tgt_pos = tgt_pos     # Positional encoding for the target
        self.projection_layer = projection_layer

    def encode(self, src, src_mask):
        """
        Performs the encoding part of the model.

        Args:
            src (torch.Tensor): The source sequence tensor.
            src_mask: The mask for the source sequence.

        Returns:
            The encoder output tensor.
        """
        src = self.src_embed(src)
        src = self.src_pos(src)
        return self.encoder(src, src_mask)

    def decode(self, encoder_output, src_mask, tgt, tgt_mask):
        """
        Performs the decoding part of the model.

        Args:
            encoder_output (torch.Tensor): The output from the encode method.
            src_mask: The mask for the source sequence.
            tgt (torch.Tensor): The target sequence tensor.
            tgt_mask: The mask for the target sequence.

        Returns:
            The decoder output tensor.
        """
        tgt = self.tgt_embed(tgt)
        tgt = self.tgt_pos(tgt)
        return self.decoder(tgt, encoder_output, src_mask, tgt_mask)

    def project(self, x):
        """
        Applies the final projection layer.
        """
        return self.projection_layer(x)

```

1.  **`encode`**: Takes the source sentence, embeds it, adds positional encoding, and passes it through the encoder stack.
2.  **`decode`**: Takes the target sentence and the encoder's output, embeds the target, adds positional encoding, and passes it all through the decoder stack.
3.  **`project`**: Takes the final decoder output and projects it into word probabilities.

We have successfully built a Transformer from the ground up. Every piece is accounted for, from the initial embedding to the final prediction.

## Credits
This lesson was provided by Gemini 2.5 Pro (AI).
