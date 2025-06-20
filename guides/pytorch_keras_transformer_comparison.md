<a name="top"></a>
# PyTorch vs. Keras/TensorFlow for Transformer Implementation

## Introduction
Transformers are a foundational architecture in modern machine learning. Two of the most popular deep learning frameworks for building them are **PyTorch** and **Keras** (with a TensorFlow backend). While both can achieve similar results, they offer different developer experiences rooted in their core design philosophies.

This document provides a side-by-side comparison of implementing Transformer components and training pipelines, highlighting:
*   The trade-offs between granular control (PyTorch) and ease of use (Keras).
*   Differences in API style and abstraction levels.
*   How common components are realized in each framework.

### Key Philosophies at a Glance
*   **PyTorch:**
    *   **Imperative & Pythonic:** Feels like writing standard Python. Operations are executed as they are called.
    *   **Explicit & Granular:** Gives the developer full control over every detail, from tensor manipulation to the training loop.
    *   **Flexibility for Research:** Its transparency makes it a favorite in research for experimenting with novel architectures and complex training schemes.

*   **Keras/TensorFlow:**
    *   **Declarative & Abstracted:** Define layers of a model and then compile it. The framework handles the underlying graph execution.
    *   **High-Level & User-Friendly:** Offers pre-built, optimized layers and a simple `model.fit()` API for rapid prototyping.
    *   **Production-Ready:** Strong support for easy deployment and scalability with the TensorFlow ecosystem (e.g., TFX, TensorFlow Serving).

---
### Table of Contents
1.  [Input Embeddings & Positional Encoding](#1-input-embeddings-and-positional-encoding)
2.  [Multi-Head Attention & Core Block](#2-multi-head-attention-mha-and-core-block)
3.  [Full Model Assembly](#3-full-model-assembly)
4.  [Data Simulation & Preparation](#4-data-simulation-and-preparation)
5.  [Training Loop](#5-training-loop)
6.  [Prediction / Inference](#6-prediction--inference)
7.  [Conclusion: Which Should You Choose?](#7-conclusion-which-framework-should-you-choose)

---

## 1. Input Embeddings and Positional Encoding
This first step converts input token IDs into information-rich vectors that combine semantic meaning with sequence order.

<details>
<summary><b>Click to view PyTorch Implementation</b></summary>

```python
import torch
import torch.nn as nn
import math

class InputEmbeddings(nn.Module):
    def __init__(self, d_model: int, vocab_size: int):
        super().__init__()
        self.d_model = d_model
        self.vocab_size = vocab_size
        self.embedding = nn.Embedding(vocab_size, d_model)

    def forward(self, x):
        # Scale embeddings as per the paper "Attention Is All You Need"
        return self.embedding(x) * math.sqrt(self.d_model)

class PositionalEncoding(nn.Module):
    def __init__(self, d_model: int, seq_len: int, dropout: float):
        super().__init__()
        self.d_model = d_model
        self.seq_len = seq_len
        self.dropout = nn.Dropout(dropout)

        # Create a matrix for positional encodings
        pe = torch.zeros(seq_len, d_model)
        position = torch.arange(0, seq_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * (-math.log(10000.0) / d_model))
        
        # Apply sine to even indices, cosine to odd
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0)

        # Register 'pe' as a buffer, not a trainable parameter
        self.register_buffer('pe', pe)

    def forward(self, x):
        # Add positional encoding to the input embeddings
        x = x + (self.pe[:, :x.shape[1], :]).requires_grad_(False)
        return self.dropout(x)
```
</details>

<details>
<summary><b>Click to view Keras Implementation</b></summary>

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
        # One embedding layer for tokens
        self.token_emb = layers.Embedding(input_dim=vocab_size, output_dim=embed_dim)
        # Another embedding layer for positions
        self.pos_emb = layers.Embedding(input_dim=maxlen, output_dim=embed_dim)

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
```
</details>

### Comparison
*   **PyTorch Approach:**
    *   Defines two separate, explicit modules: `InputEmbeddings` and `PositionalEncoding`.
    *   Pre-computes fixed sinusoidal positional encodings as described in the original paper.
    *   Uses a non-trainable `buffer` to store the positional encoding matrix.

*   **Keras Approach:**
    *   Combines both token and positional embeddings into a single, convenient custom layer.
    *   Uses a *learned* positional embedding, which is another common and effective technique.

*   **Key Differences:**
    *   **Modularity:** PyTorch uses two distinct modules; Keras encapsulates both into one layer.
    *   **Encoding Method:** PyTorch uses fixed sinusoidal encodings; Keras uses learned positional embeddings.
    *   **API Style:** PyTorch is more explicit about the math; Keras is more declarative.

---

## 2. Multi-Head Attention (MHA) and Core Block
The core of the Transformer is its MHA mechanism, followed by a Feed-Forward Network (FFN), with residual connections and layer normalization.

<details>
<summary><b>Click to view PyTorch Implementation</b></summary>

```python
# (Includes MultiHeadAttentionBlock, FeedForwardBlock, LayerNormalization, ResidualConnection)
# ... Full PyTorch code from user prompt ...
class MultiHeadAttentionBlock(nn.Module):
    def __init__(self, d_model: int, h: int, dropout: float):
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

# ... (And other helper classes like FeedForwardBlock, LayerNormalization, etc.)
```
</details>

<details>
<summary><b>Click to view Keras Implementation</b></summary>

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
        # Built-in, optimized MHA layer
        self.att = layers.MultiHeadAttention(num_heads=num_heads, key_dim=embed_dim)
        # FFN using a Sequential model
        self.ffn = keras.Sequential(
            [layers.Dense(ff_dim, activation="relu"), layers.Dense(embed_dim),]
        )
        self.layernorm1 = layers.LayerNormalization(epsilon=1e-6)
        self.layernorm2 = layers.LayerNormalization(epsilon=1e-6)
        self.dropout1 = layers.Dropout(rate)
        self.dropout2 = layers.Dropout(rate)

    def call(self, inputs, training=False):
        attn_output = self.att(inputs, inputs)
        attn_output = self.dropout1(attn_output, training=training)
        out1 = self.layernorm1(inputs + attn_output)
        ffn_output = self.ffn(out1)
        ffn_output = self.dropout2(ffn_output, training=training)
        return self.layernorm2(out1 + ffn_output)
    
    # ... get_config and from_config methods ...
```
</details>

### Comparison
*   **PyTorch Approach:**
    *   Builds MHA, FFN, and Add & Norm as separate, composable modules from scratch.
    *   Manually implements the scaled dot-product attention logic, including reshaping tensors for multiple heads.
    *   Demonstrates a Pre-LN (Layer-Normalization) style (`x + Sublayer(Norm(x))`), which is often more stable for deep Transformers.

*   **Keras Approach:**
    *   Encapsulates the entire encoder block (MHA, FFN, Add & Norm) into a single `TransformerBlock` layer.
    *   Uses the highly optimized, built-in `layers.MultiHeadAttention`, abstracting away the internal mechanics.
    *   Demonstrates a Post-LN style (`Norm(x + Sublayer(x))`).

*   **Key Differences:**
    *   **Abstraction:** PyTorch shows a detailed, low-level implementation. Keras uses a high-level, convenient built-in layer.
    *   **Integration:** PyTorch defines separate components. Keras integrates them into one larger, reusable block.
    *   **Customizability:** The PyTorch approach offers maximum flexibility for modifying the attention mechanism itself.

---

## 3. Full Model Assembly
Assembling the full model involves stacking the core blocks and adding the final projection layer.

<details>
<summary><b>Click to view PyTorch Implementation</b></summary>

```python
# (Includes EncoderBlock, DecoderBlock, Encoder, Decoder, ProjectionLayer, Transformer)
# ... Full PyTorch code from user prompt ...
class Encoder(nn.Module):
    # ...
class Decoder(nn.Module):
    # ...
class Transformer(nn.Module):
    def __init__(self, encoder: Encoder, decoder: Decoder, src_embed, tgt_embed, src_pos, tgt_pos, projection_layer):
        super().__init__()
        self.encoder = encoder
        self.decoder = decoder
        self.src_embed = src_embed
        self.tgt_embed = tgt_embed
        self.src_pos = src_pos
        self.tgt_pos = tgt_pos
        self.projection_layer = projection_layer

    def encode(self, src, src_mask):
        # ...
    def decode(self, encoder_output, src_mask, tgt, tgt_mask):
        # ...
    def project(self, x):
        # ...
```
</details>

<details>
<summary><b>Click to view Keras Implementation</b></summary>

```python
def build_model(maxlen, vocab_size, embed_dim, num_heads, ff_dim, num_transformer_blocks, dropout_rate=0.1):
    """
    Builds an encoder-only Transformer model.
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
```
</details>

### Comparison
*   **PyTorch Approach:**
    *   Builds a canonical Encoder-Decoder Transformer suitable for sequence-to-sequence tasks like translation.
    *   Highly modular, with separate classes for `EncoderBlock`, `DecoderBlock`, `Encoder`, `Decoder`, and the final `Transformer`.

*   **Keras Approach:**
    *   The example uses a `build_model` function to construct an **encoder-only** Transformer, common for tasks like sequence classification or labeling.
    *   Assembles the model using fewer, larger building blocks (`TokenAndPositionEmbedding`, `TransformerBlock`).

*   **Key Differences:**
    *   **Architecture:** The PyTorch example is a full Encoder-Decoder. The Keras example is Encoder-only. A full Encoder-Decoder could be built in Keras, but it would require more structure.
    *   **Model Definition:** PyTorch defines the full architecture in a final `Transformer` `nn.Module`. Keras uses a builder function that returns a compiled `keras.Model`.

---

## 4. Data Simulation and Preparation
A crucial step providing the raw material for training and evaluation.

<details>
<summary><b>Click to view Keras Data Prep (Conceptual)</b></summary>

```python
import tensorflow as tf
from tensorflow.keras.preprocessing.sequence import pad_sequences
import pandas as pd
import numpy as np

# This conceptual code shows the data preparation steps.
# The `simulate_data_v2.py` script would generate the input CSV.

# vocab_map = {'A': 0, 'T': 1, 'G': 2, 'C': 3, ...}

def tokenize_and_pad(sequences, vocab_map, maxlen):
    tokenized = [[vocab_map[token] for token in seq] for seq in sequences]
    padded = pad_sequences(tokenized, maxlen=maxlen, padding='post')
    return padded

# df = pd.read_csv("simulated_data.csv")
# X_train_seq = df["ancestor"].tolist()
# y_train_seq = df["descendant"].tolist()
# X_train = tokenize_and_pad(X_train_seq, vocab_map, maxlen)
# y_train = tokenize_and_pad(y_train_seq, vocab_map, maxlen)
# y_train_cat = tf.keras.utils.to_categorical(y_train, num_classes=len(vocab_map))
```
</details>

### PyTorch Perspective
The didactic PyTorch code focused on architecture, but a typical workflow would use:
*   A custom `torch.utils.data.Dataset` class to handle reading files, tokenizing, and numericalizing data.
*   A `torch.utils.data.DataLoader` to wrap the `Dataset`, providing iterable batches, shuffling, and efficient data loading.

### Comparison
*   **Completeness:** The Keras example shows a more complete data pipeline from simulation to padded tensors.
*   **Tools:** Keras uses utilities like `pad_sequences`. PyTorch's ecosystem is built around the `Dataset` and `DataLoader` classes, offering immense flexibility for custom data handling.

---

## 5. Training Loop
This is where the model learns from the data. The frameworks have starkly different approaches.

<details>
<summary><b>Click to view PyTorch Training Loop (Conceptual)</b></summary>

```python
# Conceptual PyTorch Training Loop
# model = YourTransformerModel(...)
# criterion = nn.CrossEntropyLoss()
# optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
# train_dataloader = YourDataLoader(...)
# device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
# model.to(device)

model.train() # Set model to training mode
for epoch in range(num_epochs):
    for batch in train_dataloader:
        # Move data to device
        # ...
        
        # Zero gradients
        optimizer.zero_grad()
        
        # Forward pass
        output_logits = model(...) # Actual model call
        
        # Calculate loss
        loss = criterion(output_logits.view(-1, VOCAB_SIZE), target.view(-1))
        
        # Backward pass
        loss.backward()
        
        # Update weights
        optimizer.step()
```
</details>

<details>
<summary><b>Click to view Keras Training Loop</b></summary>

```python
import tensorflow as tf
from tensorflow import keras

# model = build_model(...) # As defined previously
# X_train, y_train, X_val, y_val are prepared numpy arrays

model.compile(
    optimizer="adam", 
    loss="categorical_crossentropy", 
    metrics=["accuracy", tf.keras.metrics.Precision(), tf.keras.metrics.Recall()]
)

history = model.fit(
    X_train, y_train, 
    batch_size=64, 
    epochs=5, 
    validation_data=(X_val, y_val)
)
```
</details>

### Comparison
*   **PyTorch Approach:**
    *   Requires a **manually written training loop**.
    *   Offers complete, fine-grained control over every step: data movement, gradient zeroing, forward/backward passes, and weight updates.
    *   Maximum flexibility for research, custom logging, or complex training procedures like gradient accumulation.

*   **Keras Approach:**
    *   Abstracts the loop away into the high-level `model.fit()` method.
    *   `model.compile()` is used to configure the optimizer, loss, and metrics.
    *   Significantly more concise and easier for standard training scenarios.
    *   Behavior can be customized with `tf.keras.callbacks` (e.g., for early stopping or learning rate schedules).

---

## 6. Prediction / Inference
Using the trained model to make predictions on new data.

<details>
<summary><b>Click to view PyTorch Inference (Conceptual)</b></summary>

```python
# Conceptual PyTorch Inference
# model.load_state_dict(torch.load('model.pth'))
# model.eval() # CRITICAL: Set model to evaluation mode

# input_tensor = preprocess_function(...)

with torch.no_grad(): # Disable gradient calculations
    output_logits = model(input_tensor)

# predicted_ids = torch.argmax(output_logits, dim=-1)
# predicted_sequence = postprocess_function(predicted_ids)
```
</details>

<details>
<summary><b>Click to view Keras Inference (Conceptual)</b></summary>

```python
# with keras.utils.custom_object_scope(...):
#     loaded_model = keras.models.load_model('model.keras')

# input_data = preprocess_function(...)
# prediction_probs = loaded_model.predict(input_data)
# predicted_ids = np.argmax(prediction_probs, axis=-1)
# predicted_sequence = postprocess_function(predicted_ids)
```
</details>

### Comparison
*   **PyTorch Approach:**
    *   Requires explicit `model.eval()` to turn off layers like Dropout.
    *   Requires wrapping the forward pass in a `with torch.no_grad():` block to save memory and computation.
    *   Inference is a manual forward pass: `model(...)`.

*   **Keras Approach:**
    *   Uses the convenient `model.predict()` method, which handles setting the evaluation mode and disabling gradients automatically.
    *   Requires using a `custom_object_scope` if the model was saved with custom layers.

---

## 7. Conclusion: Which Framework Should You Choose?

There is no single "best" framework; the choice depends on your goals.

| Aspect                | PyTorch                                        | Keras / TensorFlow                               |
| --------------------- | ---------------------------------------------- | ------------------------------------------------ |
| **API Style**         | Imperative, Pythonic, explicit                 | Declarative, high-level, abstract                |
| **Control vs. Ease**  | **Maximum Control:** You write everything.     | **Maximum Ease:** Pre-built, optimized components. |
| **Best for Research** | **Excellent.** Ideal for novel architectures.  | Good, but less flexible for deep modifications.   |
| **Best for Production** | Good, with tools like TorchServe.              | **Excellent.** Mature, scalable ecosystem (TFX).   |
| **Learning Curve**    | Steeper, requires understanding DL concepts.   | Gentler, great for beginners and rapid prototyping.|
| **Debugging**         | Easier, as operations are immediate. Use `pdb`.| Can be harder due to the deferred execution graph. |

In summary, choose **PyTorch** for research, maximum flexibility, and when you want to understand and control every detail of your model. Choose **Keras/TensorFlow** for rapid prototyping, building standard architectures quickly, and when you need a clear and robust path to production deployment.

## Credits
This comparison uses PyTorch code from the "Transformer Core Concepts" document (originally provided by Gemini Pro AI) and Keras code examples that were a product of collaboration with Gemini Pro as the major author.
