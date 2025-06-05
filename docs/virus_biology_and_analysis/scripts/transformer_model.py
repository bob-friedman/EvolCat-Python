import tensorflow as tf
from tensorflow.keras import layers
import numpy as np
import math # For positional encoding

# --- Masking Helper Functions ---
def create_padding_mask(seq):
    """
    Creates a padding mask for a sequence.
    Masks positions where seq == 0 (padding token).
    Output shape: (batch_size, 1, 1, seq_len) for broadcasting with attention logits.
    """
    seq = tf.cast(tf.math.equal(seq, 0), tf.float32)
    return seq[:, tf.newaxis, tf.newaxis, :]  # (batch_size, 1, 1, seq_len)

def create_look_ahead_mask(size):
    """
    Creates a look-ahead mask for the decoder's self-attention.
    Masks future tokens.
    Output shape: (seq_len, seq_len)
    """
    mask = 1 - tf.linalg.band_part(tf.ones((size, size)), -1, 0)
    return mask  # (seq_len, seq_len)

# --- 1. Positional Encoding ---
class PositionalEncoding(layers.Layer):
    def __init__(self, position, d_model):
        super(PositionalEncoding, self).__init__()
        self.position = position
        self.d_model = d_model

        angle_rads = self.get_angles(
            np.arange(position)[:, np.newaxis],
            np.arange(d_model)[np.newaxis, :],
            d_model
        )
        angle_rads[:, 0::2] = np.sin(angle_rads[:, 0::2])
        angle_rads[:, 1::2] = np.cos(angle_rads[:, 1::2])

        self.pos_encoding = angle_rads[np.newaxis, ...]
        self.pos_encoding = tf.cast(self.pos_encoding, dtype=tf.float32)

    def get_angles(self, position, i, d_model):
        angles = 1 / np.power(10000, (2 * (i // 2)) / np.float32(d_model))
        return position * angles

    def call(self, x):
        seq_len = tf.shape(x)[1]
        return x + self.pos_encoding[:, :seq_len, :]

# --- 2. Scaled Dot-Product Attention (Standalone Function) ---
def scaled_dot_product_attention(q, k, v, mask):
    matmul_qk = tf.matmul(q, k, transpose_b=True)
    dk = tf.cast(tf.shape(k)[-1], tf.float32)
    scaled_attention_logits = matmul_qk / tf.math.sqrt(dk)

    if mask is not None:
        scaled_attention_logits += (mask * -1e9)

    attention_weights = tf.nn.softmax(scaled_attention_logits, axis=-1)
    output = tf.matmul(attention_weights, v)
    return output, attention_weights

# --- 3. Multi-Head Attention ---
class MultiHeadAttention(layers.Layer):
    def __init__(self, d_model, num_heads):
        super(MultiHeadAttention, self).__init__()
        self.num_heads = num_heads
        self.d_model = d_model
        assert d_model % self.num_heads == 0
        self.depth = d_model // self.num_heads
        self.wq = layers.Dense(d_model)
        self.wk = layers.Dense(d_model)
        self.wv = layers.Dense(d_model)
        self.dense = layers.Dense(d_model)

    def split_heads(self, x, batch_size):
        x = tf.reshape(x, (batch_size, -1, self.num_heads, self.depth))
        return tf.transpose(x, perm=[0, 2, 1, 3])

    def call(self, v, k, q, mask):
        batch_size = tf.shape(q)[0]
        q_proj = self.wq(q)
        k_proj = self.wk(k)
        v_proj = self.wv(v)

        q_split = self.split_heads(q_proj, batch_size)
        k_split = self.split_heads(k_proj, batch_size)
        v_split = self.split_heads(v_proj, batch_size)

        scaled_attention, attention_weights = scaled_dot_product_attention(q_split, k_split, v_split, mask)
        scaled_attention = tf.transpose(scaled_attention, perm=[0, 2, 1, 3])
        concat_attention = tf.reshape(scaled_attention, (batch_size, -1, self.d_model))
        output = self.dense(concat_attention)
        return output, attention_weights

# --- Helper: Point-wise Feed-Forward Network ---
def point_wise_feed_forward_network(d_model, dff):
    return tf.keras.Sequential([
        layers.Dense(dff, activation='relu'),
        layers.Dense(d_model)
    ])

# --- 4. Encoder Layer ---
class EncoderLayer(layers.Layer):
    def __init__(self, d_model, num_heads, dff, rate=0.1):
        super(EncoderLayer, self).__init__()
        self.mha = MultiHeadAttention(d_model, num_heads)
        self.ffn = point_wise_feed_forward_network(d_model, dff)
        self.layernorm1 = layers.LayerNormalization(epsilon=1e-6)
        self.layernorm2 = layers.LayerNormalization(epsilon=1e-6)
        self.dropout1 = layers.Dropout(rate)
        self.dropout2 = layers.Dropout(rate)

    def call(self, x, training, mask):
        attn_output, _ = self.mha(x, x, x, mask)
        attn_output = self.dropout1(attn_output, training=training)
        out1 = self.layernorm1(x + attn_output)
        ffn_output = self.ffn(out1)
        ffn_output = self.dropout2(ffn_output, training=training)
        out2 = self.layernorm2(out1 + ffn_output)
        return out2

# --- 5. Decoder Layer ---
class DecoderLayer(layers.Layer):
    def __init__(self, d_model, num_heads, dff, rate=0.1):
        super(DecoderLayer, self).__init__()
        self.mha1 = MultiHeadAttention(d_model, num_heads)
        self.mha2 = MultiHeadAttention(d_model, num_heads)
        self.ffn = point_wise_feed_forward_network(d_model, dff)
        self.layernorm1 = layers.LayerNormalization(epsilon=1e-6)
        self.layernorm2 = layers.LayerNormalization(epsilon=1e-6)
        self.layernorm3 = layers.LayerNormalization(epsilon=1e-6)
        self.dropout1 = layers.Dropout(rate)
        self.dropout2 = layers.Dropout(rate)
        self.dropout3 = layers.Dropout(rate)

    def call(self, x, enc_output, training, look_ahead_mask, padding_mask):
        attn1, attn_weights_block1 = self.mha1(x, x, x, look_ahead_mask)
        attn1 = self.dropout1(attn1, training=training)
        out1 = self.layernorm1(attn1 + x)

        attn2, attn_weights_block2 = self.mha2(enc_output, enc_output, out1, padding_mask)
        attn2 = self.dropout2(attn2, training=training)
        out2 = self.layernorm2(out1 + attn2)

        ffn_output = self.ffn(out2)
        ffn_output = self.dropout3(ffn_output, training=training)
        out3 = self.layernorm3(out2 + ffn_output)
        return out3, attn_weights_block1, attn_weights_block2

# --- 6. Encoder ---
class Encoder(layers.Layer):
    def __init__(self, num_layers, d_model, num_heads, dff, input_vocab_size,
                 maximum_position_encoding, rate=0.1):
        super(Encoder, self).__init__()
        self.d_model = d_model
        self.num_layers = num_layers
        self.embedding = layers.Embedding(input_vocab_size, d_model)
        self.pos_encoding = PositionalEncoding(maximum_position_encoding, self.d_model)
        self.enc_layers = [EncoderLayer(d_model, num_heads, dff, rate)
                           for _ in range(num_layers)]
        self.dropout = layers.Dropout(rate)

    def call(self, x, training, mask):
        # x shape: (batch_size, input_seq_len)
        seq_len = tf.shape(x)[1]
        x = self.embedding(x)  # (batch_size, input_seq_len, d_model)
        x *= tf.math.sqrt(tf.cast(self.d_model, tf.float32))
        x = self.pos_encoding(x)
        x = self.dropout(x, training=training)
        for i in range(self.num_layers):
            x = self.enc_layers[i](x, training, mask)
        return x  # (batch_size, input_seq_len, d_model)

# --- 7. Decoder ---
class Decoder(layers.Layer):
    def __init__(self, num_layers, d_model, num_heads, dff, target_vocab_size,
                 maximum_position_encoding, rate=0.1):
        super(Decoder, self).__init__()
        self.d_model = d_model
        self.num_layers = num_layers
        self.embedding = layers.Embedding(target_vocab_size, d_model)
        self.pos_encoding = PositionalEncoding(maximum_position_encoding, d_model)
        self.dec_layers = [DecoderLayer(d_model, num_heads, dff, rate)
                           for _ in range(num_layers)]
        self.dropout = layers.Dropout(rate)

    def call(self, x, enc_output, training, look_ahead_mask, padding_mask):
        # x shape: (batch_size, target_seq_len)
        seq_len = tf.shape(x)[1]
        attention_weights = {}
        x = self.embedding(x)  # (batch_size, target_seq_len, d_model)
        x *= tf.math.sqrt(tf.cast(self.d_model, tf.float32))
        x = self.pos_encoding(x)
        x = self.dropout(x, training=training)
        for i, layer in enumerate(self.dec_layers):
            x, block1, block2 = layer(x, enc_output, training, look_ahead_mask, padding_mask)
            attention_weights[f'decoder_layer{i+1}_block1'] = block1
            attention_weights[f'decoder_layer{i+1}_block2'] = block2
        return x, attention_weights  # x shape: (batch_size, target_seq_len, d_model)

# --- 8. Transformer ---
class Transformer(tf.keras.Model):
    def __init__(self, num_layers, d_model, num_heads, dff, input_vocab_size,
                 target_vocab_size, pe_input, pe_target, rate=0.1):
        super(Transformer, self).__init__()
        # Store parameters if needed for other methods or config, e.g. self.num_layers = num_layers
        self.encoder = Encoder(num_layers, d_model, num_heads, dff,
                               input_vocab_size, pe_input, rate)
        self.decoder = Decoder(num_layers, d_model, num_heads, dff,
                               target_vocab_size, pe_target, rate)
        self.final_layer = layers.Dense(target_vocab_size)

    def create_masks(self, inp, tar):
        # Encoder padding mask
        enc_padding_mask = create_padding_mask(inp)

        # Used in the 2nd attention block in the decoder.
        # This padding mask is used to mask the encoder outputs.
        dec_padding_mask = create_padding_mask(inp) # Mask based on input sequence padding

        # Used in the 1st attention block in the decoder.
        # It is used to pad and mask future tokens in the input received by
        # the decoder.
        look_ahead_mask = create_look_ahead_mask(tf.shape(tar)[1])
        dec_target_padding_mask = create_padding_mask(tar) # Mask based on target sequence padding
        combined_mask = tf.maximum(dec_target_padding_mask, look_ahead_mask)

        return enc_padding_mask, combined_mask, dec_padding_mask

    def call(self, inputs, training=False): # Default training to False for inference
        inp, tar = inputs # Expect a tuple (inp, tar)

        enc_padding_mask, look_ahead_combined_mask, dec_padding_mask = self.create_masks(inp, tar)

        enc_output = self.encoder(inp, training, enc_padding_mask)  # (batch_size, inp_seq_len, d_model)

        # dec_output.shape == (batch_size, tar_seq_len, d_model)
        dec_output, attention_weights = self.decoder(
            tar, enc_output, training, look_ahead_combined_mask, dec_padding_mask)

        final_output = self.final_layer(dec_output)  # (batch_size, tar_seq_len, target_vocab_size)

        return final_output, attention_weights

if __name__ == "__main__":
    print("Transformer model components defined.")

    # Parameters for tests
    batch_size_test = 1
    max_len_test = 60
    d_model_test = 128
    num_heads_test = 8
    dff_test = 512

    print(f"\n--- Testing PositionalEncoding ---")
    # ... (previous tests for PE, MHA, EncoderLayer, DecoderLayer remain unchanged) ...
    try:
        pos_enc_layer = PositionalEncoding(position=max_len_test, d_model=d_model_test)
        dummy_input_emb = tf.random.uniform((batch_size_test, max_len_test, d_model_test))
        pos_encoded_output = pos_enc_layer(dummy_input_emb)
        print(f"Positional Encoding Input Shape: {dummy_input_emb.shape}")
        print(f"Positional Encoding Output Shape: {pos_encoded_output.shape}")
        if pos_encoded_output.shape == dummy_input_emb.shape:
            print("PositionalEncoding basic test: SUCCESS")
        else:
            print("PositionalEncoding basic test: FAILED (shape mismatch)")
    except Exception as e:
        print(f"Error during PositionalEncoding test: {e}")

    print(f"\n--- Testing MultiHeadAttention & scaled_dot_product_attention ---")
    try:
        mha_layer = MultiHeadAttention(d_model=d_model_test, num_heads=num_heads_test)
        dummy_q = tf.random.uniform((batch_size_test, max_len_test, d_model_test))
        dummy_k = tf.random.uniform((batch_size_test, max_len_test, d_model_test))
        dummy_v = tf.random.uniform((batch_size_test, max_len_test, d_model_test))
        mha_output, mha_attn_weights = mha_layer(dummy_v, dummy_k, dummy_q, mask=None)
        print(f"MHA Input Q Shape: {dummy_q.shape}")
        print(f"MHA Output Shape: {mha_output.shape}")
        print(f"MHA Attention Weights Shape: {mha_attn_weights.shape}")
        if mha_output.shape == (batch_size_test, max_len_test, d_model_test) and \
           mha_attn_weights.shape == (batch_size_test, num_heads_test, max_len_test, max_len_test):
            print("MultiHeadAttention basic test: SUCCESS")
        else:
            print("MultiHeadAttention basic test: FAILED (shape mismatch)")
    except Exception as e:
        print(f"Error during MultiHeadAttention test: {e}")

    print(f"\n--- Testing EncoderLayer ---")
    encoder_output_for_decoder_test = None # Initialize for later use
    try:
        encoder_layer = EncoderLayer(d_model_test, num_heads_test, dff_test)
        dummy_input_encoder = tf.random.uniform((batch_size_test, max_len_test, d_model_test))
        encoder_output_for_decoder_test = encoder_layer(dummy_input_encoder, training=False, mask=None)
        print(f"EncoderLayer Input Shape: {dummy_input_encoder.shape}")
        print(f"EncoderLayer Output Shape: {encoder_output_for_decoder_test.shape}")
        if encoder_output_for_decoder_test.shape == dummy_input_encoder.shape:
             print("EncoderLayer basic test: SUCCESS")
        else:
            print("EncoderLayer basic test: FAILED (shape mismatch)")
    except Exception as e:
        print(f"Error during EncoderLayer test: {e}")

    print(f"\n--- Testing DecoderLayer ---")
    try:
        decoder_layer = DecoderLayer(d_model_test, num_heads_test, dff_test)
        dummy_input_decoder = tf.random.uniform((batch_size_test, max_len_test, d_model_test))
        enc_output_dummy = encoder_output_for_decoder_test if encoder_output_for_decoder_test is not None else tf.random.uniform((batch_size_test, max_len_test, d_model_test))
        decoder_output, dec_attn1, dec_attn2 = decoder_layer(
            dummy_input_decoder, enc_output_dummy, training=False,
            look_ahead_mask=None, padding_mask=None)
        print(f"DecoderLayer Input Shape: {dummy_input_decoder.shape}")
        print(f"DecoderLayer Encoder Output Shape (for MHA2): {enc_output_dummy.shape}")
        print(f"DecoderLayer Output Shape: {decoder_output.shape}")
        print(f"DecoderLayer Attn1 (Self-Attn) Shape: {dec_attn1.shape}")
        print(f"DecoderLayer Attn2 (Enc-Dec Attn) Shape: {dec_attn2.shape}")
        if decoder_output.shape == dummy_input_decoder.shape and \
           dec_attn1.shape == (batch_size_test, num_heads_test, max_len_test, max_len_test) and \
           dec_attn2.shape == (batch_size_test, num_heads_test, max_len_test, max_len_test):
            print("DecoderLayer basic test: SUCCESS")
        else:
            print("DecoderLayer basic test: FAILED (shape or attention weights shape mismatch)")
    except Exception as e:
        print(f"Error during DecoderLayer test: {e}")

    print(f"\n--- Testing Transformer (Full Model) ---")
    try:
        num_layers_transformer = 2
        input_vocab_size_transformer = 7
        target_vocab_size_transformer = 7
        max_pos_enc_input = max_len_test
        max_pos_enc_target = max_len_test

        transformer_model = Transformer(
            num_layers_transformer, d_model_test, num_heads_test, dff_test,
            input_vocab_size_transformer, target_vocab_size_transformer,
            pe_input=max_pos_enc_input, pe_target=max_pos_enc_target
        )

        dummy_input_seq = tf.random.uniform((batch_size_test, max_len_test), dtype=tf.int64, minval=0, maxval=input_vocab_size_transformer-1)
        dummy_target_seq = tf.random.uniform((batch_size_test, max_len_test), dtype=tf.int64, minval=0, maxval=target_vocab_size_transformer-1)

        predictions, attn_weights = transformer_model((dummy_input_seq, dummy_target_seq), training=False)
        print(f"Transformer Input Seq Shape: {dummy_input_seq.shape}")
        print(f"Transformer Target Seq Shape: {dummy_target_seq.shape}")
        print(f"Transformer Predictions Shape: {predictions.shape}")
        print(f"Transformer Attention Weights Keys: {list(attn_weights.keys())}")

        expected_pred_shape = (batch_size_test, max_len_test, target_vocab_size_transformer)
        if predictions.shape == expected_pred_shape:
            print("Transformer basic test: SUCCESS (shapes match)")
        else:
            print(f"Transformer basic test: FAILED (shape mismatch). Expected {expected_pred_shape}, Got {predictions.shape}")

    except Exception as e:
        print(f"Error during Transformer test: {e}")
```
