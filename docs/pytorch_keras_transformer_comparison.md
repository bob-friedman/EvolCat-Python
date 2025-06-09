# PyTorch vs. Keras/TensorFlow for Transformer Implementation

## Introduction
Transformers have become a foundational architecture in modern machine learning, powering advancements in natural language processing, computer vision, and more. Two of the most popular deep learning frameworks for building Transformers are PyTorch and Keras (often with a TensorFlow backend). While both can achieve similar results, they offer different developer experiences due to their design philosophies.

This document provides a side-by-side comparison of implementing core Transformer components and training pipelines using PyTorch and Keras. The PyTorch examples are drawn from the didactic `transformer_core_concepts.md` guide, emphasizing a from-scratch understanding. The Keras examples utilize user-provided code snippets that demonstrate a more integrated and application-focused approach.

The goal is to highlight:
- Differences in API style and abstraction levels.
- How common Transformer components are realized in each framework.
- The trade-offs between granular control (PyTorch) and ease of use/rapid prototyping (Keras).

Understanding these distinctions can help developers choose the framework that best suits their project needs, whether it's deep research into model internals or building and deploying robust ML applications efficiently.

## 1. Input Embeddings and Positional Encoding
### PyTorch Implementation
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
### Keras Implementation
```python
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers

class TokenAndPositionEmbedding(layers.Layer):
    def __init__(self, maxlen, vocab_size, embed_dim, **kwargs):
        super().__init__(**kwargs)
        self.maxlen = maxlen
        self.vocab_size = vocab_size
        self.embed_dim = embed_dim
        self.token_emb = layers.Embedding(input_dim=vocab_size, output_dim=embed_dim)
        self.pos_emb = layers.Embedding(input_dim=maxlen, output_dim=embed_dim)

    def build(self, input_shape):
        # Potentially add weights or other initializations here if needed
        # For example, self.add_weight(...)
        super().build(input_shape) # Be sure to call this at the end

    def call(self, x):
        maxlen = tf.shape(x)[-1]
        positions = tf.range(start=0, limit=maxlen, delta=1)
        positions = self.pos_emb(positions)
        x = self.token_emb(x)
        return x + positions

    def get_config(self):
        config = super().get_config()
        config.update({
            "maxlen": self.maxlen,
            "vocab_size": self.vocab_size,
            "embed_dim": self.embed_dim,
        })
        return config

    @classmethod
    def from_config(cls, config):
        return cls(**config)
```
### Comparison
Both PyTorch and Keras approaches achieve the same goal: converting input token IDs into information-rich vectors that combine semantic meaning (from token embeddings) with sequence order information (from positional embeddings).

**PyTorch Approach:**
- Explicitly defines two separate modules: `InputEmbeddings` for token embedding and `PositionalEncoding` for adding positional information.
- The `InputEmbeddings` module uses `nn.Embedding` and scales the embeddings by `sqrt(d_model)`.
- The `PositionalEncoding` module pre-computes sinusoidal encoding values and adds them to the token embeddings. It uses a registered buffer (`self.register_buffer`) for the positional encoding matrix, as these are not learned parameters.
- The two steps are applied sequentially to the input tensor.

**Keras Approach:**
- Combines token and positional embedding into a single custom layer: `TokenAndPositionEmbedding(layers.Layer)`.
- It initializes two `layers.Embedding` sub-layers internally: one for tokens (`self.token_emb`) and one for positions (`self.pos_emb`).
- In the `call` method, it generates position indices using `tf.range` and then embeds both the input tokens and these position indices.
- The token embeddings and position embeddings are then added together.
- Keras layers are encouraged to implement `build()` and `get_config()` for proper serialization and adherence to the Keras API.

**Key Differences:**
- **Modularity:** PyTorch uses two distinct, sequential modules. Keras encapsulates both functionalities within a single, reusable layer.
- **Positional Encoding Method:** The PyTorch example uses fixed sinusoidal positional encodings as described in the original Transformer paper. The Keras example uses a *learned* positional embedding layer (`self.pos_emb`), which is another common technique.
- **API Style:** PyTorch code tends to be more explicit about the mathematical operations and module composition. Keras abstracts some of this into its Layer class, promoting a more declarative style.

## 2. Multi-Head Attention (MHA) and Core Block
### PyTorch Implementation
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
### Keras Implementation
```python
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers

class TransformerBlock(layers.Layer):
    def __init__(self, embed_dim, num_heads, ff_dim, rate=0.1, **kwargs):
        super().__init__(**kwargs)
        self.embed_dim = embed_dim
        self.num_heads = num_heads
        self.ff_dim = ff_dim
        self.rate = rate
        self.att = layers.MultiHeadAttention(num_heads=num_heads, key_dim=embed_dim)
        self.ffn = keras.Sequential(
            [layers.Dense(ff_dim, activation="relu"), layers.Dense(embed_dim),]
        )
        self.layernorm1 = layers.LayerNormalization(epsilon=1e-6)
        self.layernorm2 = layers.LayerNormalization(epsilon=1e-6)
        self.dropout1 = layers.Dropout(rate)
        self.dropout2 = layers.Dropout(rate)

    def call(self, inputs, training=False): # Added training argument
        attn_output = self.att(inputs, inputs)
        attn_output = self.dropout1(attn_output, training=training) # Pass training
        out1 = self.layernorm1(inputs + attn_output)
        ffn_output = self.ffn(out1)
        ffn_output = self.dropout2(ffn_output, training=training) # Pass training
        return self.layernorm2(out1 + ffn_output)

    def get_config(self):
        config = super().get_config()
        config.update({
            "embed_dim": self.embed_dim,
            "num_heads": self.num_heads,
            "ff_dim": self.ff_dim,
            "rate": self.rate,
        })
        return config

    @classmethod
    def from_config(cls, config):
        return cls(**config)
```
### Comparison
The core of the Transformer architecture lies in its Multi-Head Attention (MHA) mechanism, typically followed by a Feed-Forward Network (FFN), with residual connections and layer normalization applied around each.

**PyTorch Approach:**
- The PyTorch example builds the MHA block (`MultiHeadAttentionBlock`) from scratch. This involves:
    - Linear layers (`nn.Linear`) for creating Query (Q), Key (K), and Value (V) projections, and for the final output projection.
    - Manually reshaping Q, K, V tensors to split them across multiple attention heads.
    - A static `attention` method that computes scaled dot-product attention scores (`(Q @ K.T) / sqrt(d_k)`), applies masking, and then applies softmax.
    - Concatenating the outputs of the heads.
- A separate `FeedForwardBlock` is defined, usually consisting of two linear layers with a ReLU activation in between.
- Custom `LayerNormalization` and `ResidualConnection` modules are implemented to apply `Add & Norm` steps. The `ResidualConnection` module in the PyTorch example uses a Pre-LN style normalization (Norm -> Sublayer -> Dropout -> Add).
- These components would be explicitly combined in an `EncoderBlock` or `DecoderBlock` (as seen in `transformer_core_concepts.md`).

**Keras Approach:**
- The Keras example uses a single `TransformerBlock(layers.Layer)` which encapsulates MHA, FFN, Layer Normalization, and Dropout.
- It leverages the built-in `layers.MultiHeadAttention` layer. This Keras layer internally handles the creation of Q, K, V projections, splitting/combining heads, and the attention calculation. You primarily configure it with `num_heads` and `key_dim`.
- The Feed-Forward Network (`self.ffn`) is defined using `keras.Sequential` with `layers.Dense`.
- It uses built-in `layers.LayerNormalization` and `layers.Dropout`.
- The `call` method of `TransformerBlock` directly applies the sequence: MHA -> Dropout -> Add & Norm (with `self.layernorm1`) -> FFN -> Dropout -> Add & Norm (with `self.layernorm2`). This is a Post-LN style for the first sublayer and Post-LN for the second as well, though the exact order of operations (norm first vs. sublayer first) can vary in Transformer implementations.

**Key Differences:**
- **Level of Abstraction for MHA:** PyTorch shows a manual, detailed implementation of MHA, offering fine-grained control and understanding of each step. Keras uses a high-level `layers.MultiHeadAttention` layer, abstracting away the internal mechanics.
- **Component Integration:** The PyTorch example defines MHA, FFN, and Add & Norm as separate, composable modules. The Keras `TransformerBlock` integrates these into a single, larger layer representing a standard Transformer encoder block.
- **Normalization Style:** The PyTorch example explicitly implemented Pre-LN style residual connections (`x + self.dropout(sublayer(self.norm(x)))`). The Keras `TransformerBlock` shows a common pattern of `out1 = self.layernorm1(inputs + attention_output)` which is a Post-LN style for the attention sublayer, and `out2 = self.layernorm2(out1 + ffn_output)` for the FFN sublayer. Both Pre-LN and Post-LN are valid approaches with different training dynamics.
- **Customizability vs. Convenience:** PyTorch provides maximum flexibility for custom modifications to any part of the attention mechanism. Keras offers convenience and robustness with its pre-built, optimized layers, suitable for rapid prototyping and standard implementations.

## 3. Full Model Assembly
### PyTorch Implementation
```python
import torch
import torch.nn as nn
import math

# Assuming InputEmbeddings, PositionalEncoding, MultiHeadAttentionBlock,
# FeedForwardBlock, LayerNormalization, ResidualConnection are defined as in previous sections.

class EncoderBlock(nn.Module):
    def __init__(self, self_attention_block: MultiHeadAttentionBlock, feed_forward_block: FeedForwardBlock, dropout: float):
        super().__init__()
        self.self_attention_block = self_attention_block
        self.feed_forward_block = feed_forward_block
        self.residual_connections = nn.ModuleList([ResidualConnection(dropout) for _ in range(2)])

    def forward(self, x, src_mask):
        x = self.residual_connections[0](x, lambda x_att: self.self_attention_block(x_att, x_att, x_att, src_mask))
        x = self.residual_connections[1](x, self.feed_forward_block)
        return x

class DecoderBlock(nn.Module):
    def __init__(self, self_attention_block: MultiHeadAttentionBlock, cross_attention_block: MultiHeadAttentionBlock, feed_forward_block: FeedForwardBlock, dropout: float):
        super().__init__()
        self.self_attention_block = self_attention_block
        self.cross_attention_block = cross_attention_block
        self.feed_forward_block = feed_forward_block
        self.residual_connections = nn.ModuleList([ResidualConnection(dropout) for _ in range(3)])

    def forward(self, x, encoder_output, src_mask, tgt_mask):
        x = self.residual_connections[0](x, lambda x_att: self.self_attention_block(x_att, x_att, x_att, tgt_mask))
        x = self.residual_connections[1](x, lambda x_cross: self.cross_attention_block(x_cross, encoder_output, encoder_output, src_mask))
        x = self.residual_connections[2](x, self.feed_forward_block)
        return x

class Encoder(nn.Module):
    def __init__(self, layers: nn.ModuleList):
        super().__init__()
        self.layers = layers
        self.norm = LayerNormalization() # Assuming LayerNormalization is defined

    def forward(self, x, mask):
        for layer in self.layers:
            x = layer(x, mask)
        return self.norm(x)

class Decoder(nn.Module):
    def __init__(self, layers: nn.ModuleList):
        super().__init__()
        self.layers = layers
        self.norm = LayerNormalization() # Assuming LayerNormalization is defined

    def forward(self, x, encoder_output, src_mask, tgt_mask):
        for layer in self.layers:
            x = layer(x, encoder_output, src_mask, tgt_mask)
        return self.norm(x)

class ProjectionLayer(nn.Module):
    def __init__(self, d_model: int, vocab_size: int):
        super().__init__()
        self.proj = nn.Linear(d_model, vocab_size)

    def forward(self, x):
        return torch.log_softmax(self.proj(x), dim=-1)

class Transformer(nn.Module):
    def __init__(self, encoder: Encoder, decoder: Decoder, src_embed: InputEmbeddings, tgt_embed: InputEmbeddings, src_pos: PositionalEncoding, tgt_pos: PositionalEncoding, projection_layer: ProjectionLayer):
        super().__init__()
        self.encoder = encoder
        self.decoder = decoder
        self.src_embed = src_embed
        self.tgt_embed = tgt_embed
        self.src_pos = src_pos
        self.tgt_pos = tgt_pos
        self.projection_layer = projection_layer

    def encode(self, src, src_mask):
        src = self.src_embed(src)
        src = self.src_pos(src)
        return self.encoder(src, src_mask)

    def decode(self, encoder_output, src_mask, tgt, tgt_mask):
        tgt = self.tgt_embed(tgt)
        tgt = self.tgt_pos(tgt)
        return self.decoder(tgt, encoder_output, src_mask, tgt_mask)

    def project(self, x):
        return self.projection_layer(x)
```
### Keras Implementation
```python
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers

# Assuming TokenAndPositionEmbedding and TransformerBlock are defined
# as in previous sections or available in the environment.
# For clarity, their definitions are usually in separate files
# and imported, e.g.:
# from .transformer_model import TokenAndPositionEmbedding, TransformerBlock

# class TokenAndPositionEmbedding(layers.Layer): ... (see section 1)
# class TransformerBlock(layers.Layer): ... (see section 2)

def build_model(maxlen, vocab_size, embed_dim, num_heads, ff_dim, num_transformer_blocks, dropout_rate=0.1):
    """
    Builds a Transformer-based model.
    Note: This example is more like an encoder-only Transformer for sequence labeling
          or auto-encoding rather than a full encoder-decoder for translation.
    """
    inputs = layers.Input(shape=(maxlen,))
    embedding_layer = TokenAndPositionEmbedding(maxlen, vocab_size, embed_dim)
    x = embedding_layer(inputs)

    for _ in range(num_transformer_blocks):
        transformer_block = TransformerBlock(embed_dim, num_heads, ff_dim, rate=dropout_rate)
        x = transformer_block(x)

    outputs = layers.Dense(vocab_size, activation="softmax")(x)

    model = keras.Model(inputs=inputs, outputs=outputs)
    return model

# Example usage:
# vocab_size = 20000  # Only consider the top 20k words
# maxlen = 200  # Max sequence length
# embed_dim = 256  # Embedding size for each token
# num_heads = 2  # Number of attention heads
# ff_dim = 256  # Hidden layer size in feed forward network inside transformer
# num_transformer_blocks = 1 # Number of transformer blocks

# model = build_model(
#     maxlen=maxlen,
#     vocab_size=vocab_size,
#     embed_dim=embed_dim,
#     num_heads=num_heads,
#     ff_dim=ff_dim,
#     num_transformer_blocks=num_transformer_blocks
# )
```
### Comparison
Assembling the full Transformer model involves bringing together the embedding layers, encoder stack, decoder stack (for sequence-to-sequence tasks), and a final projection layer.

**PyTorch Approach:**
- The PyTorch example defines separate classes for `EncoderBlock` and `DecoderBlock`, which are the building blocks for the encoder and decoder stacks respectively.
- The `Encoder` class takes a list of `EncoderBlock` instances and applies them sequentially, followed by a final `LayerNormalization`.
- Similarly, the `Decoder` class takes a list of `DecoderBlock` instances.
- A `ProjectionLayer` is defined, which is a simple linear layer followed by a softmax, to map the decoder's output to vocabulary probabilities.
- The final `Transformer` class then composes these: source and target embeddings, source and target positional encodings, the `Encoder`, the `Decoder`, and the `ProjectionLayer`. It defines `encode`, `decode`, and `project` methods to manage the data flow.
- This approach is highly modular, allowing easy modification or replacement of any component (e.g., using a different type of attention, adding more layers).

**Keras Approach:**
- The Keras example uses a `build_model` function to construct the model. This function is more analogous to an encoder-only Transformer, as it doesn't explicitly define a separate decoder stack for a sequence-to-sequence task (the provided Keras example is geared towards sequence-to-sequence on the same vocabulary, like `A->G` transition, effectively an auto-encoder or sequence labeling task).
- It defines an `Input` layer.
- It uses the `TokenAndPositionEmbedding` layer (defined earlier) for the input.
- It then passes the output through one or more `TransformerBlock` layers (defined earlier). The `TransformerBlock` itself represents a complete encoder layer (MHA + FFN + Add & Norm).
- A final `layers.Dense` layer with a `softmax` activation acts as the projection layer to output vocabulary probabilities.
- The entire model is wrapped in a `keras.Model` object.

**Key Differences:**
- **Modularity and Granularity:** The PyTorch example breaks the Transformer down into more numerous, finer-grained modules (`EncoderBlock`, `DecoderBlock`, `Encoder`, `Decoder`, `ProjectionLayer`, `Transformer`). Keras, using the `TransformerBlock` which is already an entire encoder layer, assembles the model with fewer, larger building blocks.
- **Encoder/Decoder Structure:** The PyTorch `Transformer` class explicitly includes both an encoder and a decoder, making it suitable for canonical sequence-to-sequence tasks (e.g., translation). The Keras `build_model` function, as provided, constructs an architecture that processes an input sequence and produces an output sequence of the same length, which is common in tasks like sequence labeling or auto-encoding, rather than a full encoder-decoder for distinct source/target sequences. To make the Keras model a true encoder-decoder, one would typically define separate `Encoder` and `Decoder` stacks using `TransformerBlock` (or similar) and then connect them, similar to the PyTorch structure.
- **Model Definition:** In PyTorch, the `Transformer` class itself is an `nn.Module` defining the full architecture. In Keras, the `build_model` function returns a compiled `keras.Model` instance.
- **Flexibility for Research:** The PyTorch structure is often favored for research due to its explicitness, making it easier to experiment with novel architectural changes within any component. Keras provides a more streamlined path for building standard architectures.

## 4. Data Simulation and Preparation
### Overview
The Keras example includes a script `simulate_data_v2.py` (provided by the user) designed to generate a toy dataset. Its purpose is to create biologically-inspired ancestor and descendant DNA-like sequences with a specific mutational pattern (e.g., 'A' -> 'G', with some sites becoming 'R' for polymorphism). This data is saved to a CSV file, which is then used by the training script. This step is crucial for any machine learning pipeline, as it provides the raw material for training and evaluation.
### Keras Implementation
```python
import tensorflow as tf
from tensorflow.keras.preprocessing.sequence import pad_sequences
import numpy as np
# import simulate_data_v2 # Assuming this module is available

# From train.py (conceptual)
# vocab_map = {nt: i for i, nt in enumerate(simulate_data_v2.FULL_VOCABULARY)}
# inv_vocab_map = {i: nt for nt, i in vocab_map.items()}

def tokenize_and_pad(sequences, vocab_map, maxlen):
    tokenized_sequences = [[vocab_map[token] for token in seq] for seq in sequences]
    padded_sequences = pad_sequences(tokenized_sequences, maxlen=maxlen, padding='post', truncating='post')
    return padded_sequences

# Example usage from train.py (conceptual)
# df = pd.read_csv("simulated_data.csv")
# X_train_seq = df["ancestor"].tolist()
# y_train_seq = df["descendant"].tolist()

# X_train = tokenize_and_pad(X_train_seq, vocab_map, maxlen)
# y_train = tokenize_and_pad(y_train_seq, vocab_map, maxlen)
# y_train = tf.keras.utils.to_categorical(y_train, num_classes=len(vocab_map)) # For categorical crossentropy
```
### PyTorch Perspective
The PyTorch Transformer code in `transformer_core_concepts.md` was primarily didactic, focusing on the model architecture itself. It did not include explicit data loading or preprocessing steps. In a typical PyTorch workflow:
- You would use or create a `Dataset` class (from `torch.utils.data.Dataset`) to load and process your data. This class would handle reading data from files (like the CSV generated by `simulate_data_v2.py`), tokenizing sequences, converting tokens to numerical IDs, and perhaps applying initial padding or truncation.
- A `DataLoader` (from `torch.utils.data.DataLoader`) would then wrap the `Dataset` to provide iterable batches of data to the model during training, handling shuffling, batching, and potentially multiprocessing for efficiency.
- Tokenization would involve creating a vocabulary mapping similar to `vocab_map` in the Keras example.
- Padding would often be handled either within the `Dataset`'s `__getitem__` or by a custom `collate_fn` passed to the `DataLoader`, ensuring all sequences in a batch have the same length, often by padding to the maximum length in that batch or a predefined `maxlen`.
### Comparison
**Keras Approach:**
- The Keras example demonstrates a more complete pipeline where data simulation (`simulate_data_v2.py`) and data preparation (tokenization, padding in `train.py`) are explicit steps.
- It uses pandas to read the CSV, then custom Python for tokenization based on a `vocab_map`, and `tf.keras.preprocessing.sequence.pad_sequences` for padding.

**PyTorch Approach (General Practice):**
- While not shown in the didactic model code, PyTorch emphasizes flexible data handling through its `Dataset` and `DataLoader` classes.
- This allows for complex data loading and preprocessing logic to be encapsulated and customized.
- Padding strategies can be very flexible (e.g., dynamic padding per batch).

**Key Differences:**
- **Completeness of Example:** The Keras example provides a runnable data pipeline. The PyTorch architectural example would require separate implementation of data loading and preprocessing.
- **Tools:** Keras often uses its preprocessing utilities like `pad_sequences`. PyTorch relies on its `Dataset`/`DataLoader` paradigm, often augmented with custom logic or libraries like `torchtext` (though `torchtext` is evolving).
- **Integration:** In Keras, data preprocessing can sometimes be integrated as Keras layers themselves (e.g., `TextVectorization` layer for tokenization), though the example here uses offline preprocessing. PyTorch's `Dataset` objects are central to its training loop.

## 5. Training Loop
### PyTorch Implementation (Conceptual)
The PyTorch `transformer_core_concepts.md` document focused on model architecture and did not include a training loop. However, a standard PyTorch training loop typically involves the following manual steps:
```python
# Conceptual PyTorch Training Loop Snippet
# import torch
# import torch.nn as nn
# model = YourTransformerModel(...) # Your model definition
# criterion = nn.CrossEntropyLoss(...) # Or another suitable loss
# optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
# num_epochs = 10
# train_dataloader = YourDataLoader(...) # Your data loader instance
# device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
# model.to(device)
# VOCAB_SIZE = model.projection_layer.proj.out_features # Example to get vocab_size

# model.train() # Set model to training mode
# for epoch in range(num_epochs):
#     for batch in train_dataloader: # Assuming train_dataloader yields batches of data
#         # Assuming batch is a tuple or object with src, tgt_input, tgt_output, src_mask, tgt_mask
#         source_seq = batch.src.to(device)
#         target_seq_input = batch.tgt_input.to(device) # e.g., <sos> + sequence
#         target_seq_output = batch.tgt_output.to(device) # e.g., sequence + <eos>
#         src_mask = batch.src_mask.to(device) # Appropriate masks
#         tgt_mask = batch.tgt_mask.to(device) # Appropriate masks
#
#         # Zero gradients
#         optimizer.zero_grad()
#
#         # Forward pass
#         # For a full Transformer (encoder-decoder):
#         # output_logits = model(source_seq, target_seq_input, src_mask, tgt_mask)
#         # For an encoder-only model (like the Keras example structure):
#         # output_logits = model.project(model.encode(source_seq, src_mask)) # If model has encode & project methods
#         # Actual call depends on your model's forward() method for the specific task
#         # For this conceptual snippet, let's assume output_logits is the direct model output
#         output_logits = model(source_seq, target_seq_input, src_mask, tgt_mask) # Placeholder for actual call
#
#         # Calculate loss (requires reshaping logits or targets for CrossEntropyLoss)
#         # output_logits shape: (batch_size, seq_len, vocab_size)
#         # target_seq_output shape: (batch_size, seq_len)
#         loss = criterion(output_logits.reshape(-1, VOCAB_SIZE), target_seq_output.reshape(-1))
#
#         # Backward pass
#         loss.backward()
#
#         # Update weights
#         optimizer.step()
#
#     # print(f"Epoch {epoch+1}/{num_epochs}, Loss: {loss.item()}")
```
This provides full control over every aspect of training, including gradient accumulation, custom learning rate schedules, and detailed logging.
### Keras Implementation
```python
import tensorflow as tf
from tensorflow import keras
# Assuming X_train, y_train, X_val, y_val are prepared
# model = build_model(...) # As defined in the Full Model Assembly section

# model.compile(optimizer="adam", loss="categorical_crossentropy", metrics=["accuracy"])

# history = model.fit(
#     X_train, y_train,
#     batch_size=64,
#     epochs=5, # Or more
#     validation_data=(X_val, y_val)
# )

# Example from user-provided train.py:
# model.compile(optimizer="adam", loss="categorical_crossentropy", metrics=["accuracy", tf.keras.metrics.Precision(), tf.keras.metrics.Recall()])
# history = model.fit(X_train, y_train, batch_size=BATCH_SIZE, epochs=EPOCHS, validation_data=(X_val, y_val))
```
### Comparison
**PyTorch Approach (Conceptual):**
- Requires manually writing the training loop.
- Explicitly define the loss function (e.g., `nn.CrossEntropyLoss`).
- Explicitly define the optimizer (e.g., `torch.optim.Adam`).
- Inside the loop, you manually:
    - Set the model to training mode (`model.train()`).
    - Iterate over data batches (typically from a `DataLoader`).
    - Move data to the appropriate device (CPU/GPU).
    - Zero out gradients (`optimizer.zero_grad()`).
    - Perform the forward pass through the model.
    - Calculate the loss.
    - Perform the backward pass to compute gradients (`loss.backward()`).
    - Update model parameters (`optimizer.step()`).
- Offers maximum flexibility for custom training procedures, complex logging, or non-standard gradient manipulations.

**Keras Approach:**
- Uses the `model.compile()` method to specify the optimizer, loss function, and metrics. Keras often infers necessary argument shapes or uses string identifiers for common optimizers/losses.
- Uses the `model.fit()` method to perform the training. This single method handles:
    - Iterating over epochs and batches.
    - Performing forward and backward passes.
    - Calculating loss and metrics.
    - Updating model weights.
    - Progress bar display and basic logging.
- Callbacks (`tf.keras.callbacks`) can be used to customize behavior during training (e.g., learning rate scheduling, early stopping, model checkpointing) without rewriting the loop.

**Key Differences:**
- **Verbosity and Control:** PyTorch is more verbose but offers fine-grained control over every step. Keras is more concise, abstracting the loop details into `fit()`.
- **Ease of Use:** Keras `fit()` is generally easier for standard training scenarios. PyTorch requires more boilerplate code but is powerful for research and complex setups.
- **Customization:** While Keras `fit()` can be customized with callbacks or by overriding `train_step`, PyTorch's manual loop is inherently more customizable.
- **Loss Functions & Optimizers:** Both offer a wide range of standard options. PyTorch often requires more explicit handling of input/output shapes for loss functions.

## 6. Prediction/Inference
### PyTorch Implementation (Conceptual)
Performing inference with a PyTorch model typically involves:
```python
# Conceptual PyTorch Inference Snippet
# import torch
# model = YourTransformerModel(...) # Your model definition
# model.load_state_dict(torch.load('model_checkpoint.pth'))
# model.eval() # Set model to evaluation mode (disables dropout, etc.)
# device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
# model.to(device)

# input_text_string = "Your input sequence here"
# # vocab_map = Your vocabulary mapping
# # id_to_char = {v: k for k, v in vocab_map.items()}
# # Preprocess input_text_string to input_tensor with tokenization, numericalization, padding
# # input_tensor = preprocess_function(input_text_string, vocab_map, maxlen).to(device)
# # input_mask = (input_tensor != padding_idx).to(device) # Create mask if needed

# with torch.no_grad(): # Disable gradient calculations for inference
#     # The actual call depends on your model's forward signature
#     # E.g., for an encoder-only model:
#     # output_logits = model.project(model.encode(input_tensor, input_mask))
#     # E.g., for a full transformer:
#     # output_logits = model(input_tensor, target_dummy_input, input_mask, target_mask)
#     output_logits = model(input_tensor) # Placeholder for actual model call

#     predicted_probabilities = torch.softmax(output_logits, dim=-1)
#     predicted_ids = torch.argmax(predicted_probabilities, dim=-1)

# # Convert predicted IDs back to characters/tokens
# # predicted_sequence = "".join([id_to_char[idx.item()] for idx in predicted_ids.squeeze()])
# # print(f"Input:      {input_text_string}")
# # print(f"Prediction: {predicted_sequence}")
```
Key steps include loading trained weights, setting the model to `eval()` mode, ensuring input tensors are on the correct device, and running the forward pass within a `torch.no_grad()` context to save memory and computation.
### Keras Implementation
```python
import tensorflow as tf
from tensorflow import keras
import numpy as np
# Assuming TokenAndPositionEmbedding and TransformerBlock are defined and model is saved as 'transformer_model.keras'
# Also assuming vocab_map and inv_vocab_map are available from training script.
# maxlen = # max sequence length used during training

# Load the model with custom objects
# with keras.utils.custom_object_scope({
#     'TokenAndPositionEmbedding': TokenAndPositionEmbedding,
#     'TransformerBlock': TransformerBlock
# }):
#     loaded_model = keras.models.load_model('transformer_model.keras')

# Example prediction function (conceptual, based on user-provided structure)
# def predict_sequence(input_text, model, vocab_map, inv_vocab_map, maxlen):
#     tokenized_input = [[vocab_map[token] for token in input_text]]
#     padded_input = tf.keras.preprocessing.sequence.pad_sequences(tokenized_input, maxlen=maxlen, padding='post', truncating='post')
#     prediction_probs = model.predict(padded_input)
#     predicted_ids = np.argmax(prediction_probs, axis=-1)

#     predicted_sequence = "".join([inv_vocab_map[idx] for idx in predicted_ids[0] if idx in inv_vocab_map])
#     return predicted_sequence

# input_dna = "ATGCGTAGCATGC" # Example input
# predicted_dna = predict_sequence(input_dna, loaded_model, vocab_map, inv_vocab_map, maxlen)
# print(f"Input:      {input_dna}")
# print(f"Prediction: {predicted_dna}")
```
### Comparison
**PyTorch Approach (Conceptual):**
- Load trained model weights using `model.load_state_dict()`.
- Set the model to evaluation mode using `model.eval()`. This is crucial as it deactivates layers like Dropout and BatchNorm, which behave differently during training and inference.
- Preprocess raw input into the required tensor format (tokenization, numericalization, padding, batching).
- Perform the forward pass within a `with torch.no_grad():` block to disable gradient computation, which is unnecessary for inference and saves resources.
- Post-process the model's output (e.g., applying softmax if logits are returned, using `argmax` to get predicted token IDs, converting IDs back to human-readable text).

**Keras Approach:**
- Load the saved model using `keras.models.load_model()`. If the model contains custom layers (like `TokenAndPositionEmbedding` or `TransformerBlock`), these must be provided via the `custom_objects` argument or a `keras.utils.custom_object_scope`.
- Preprocess raw input into the required numerical format (tokenization, padding), similar to the training data preparation.
- Use the `model.predict()` method to get output predictions. This method handles batching and the forward pass efficiently.
- Post-process the predictions (e.g., `argmax` on probabilities, decoding IDs to text).

**Key Differences:**
- **Mode Setting:** PyTorch requires explicit `model.eval()`. Keras handles this implicitly within `predict()` for its built-in layers, but for custom layers, their behavior in training vs. inference needs to be correctly defined (e.g. using the `training` argument in `call`).
- **Gradient Handling:** PyTorch needs `with torch.no_grad()`. Keras `predict()` handles this automatically.
- **Prediction Method:** Keras provides a convenient `model.predict()` method. In PyTorch, you directly call `model(...)` for the forward pass.
- **Custom Objects:** Keras requires explicit handling of custom layers during model loading. PyTorch model loading primarily deals with state dictionaries, and the class definitions must be available.

## Credits
This comparison uses PyTorch code from the "Transformer Core Concepts" document (originally provided by Gemini 2.5 Pro AI) and the Keras code examples a product of collaboration with Gemini 2.5 Pro as the major author.
