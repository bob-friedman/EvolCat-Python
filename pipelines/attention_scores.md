# Saliency Scores and Transformer Attention: A Lens for Scientific Discovery

In machine learning, a **saliency score** highlights the importance of individual input features in a model's output. It answers: "Which parts of the input did the model focus on most?" For sequences like proteins, this means identifying key amino acid residues. This creates a "saliency map," offering a transparent view into the model's decision-making.

The **attention mechanism**, core to transformer architectures, naturally provides these saliency scores. It allows a model to weigh the importance of different elements in an input sequence. For each token (e.g., an amino acid), attention calculates scores quantifying focus on other tokens. High attention between two positions suggests a strong relationship found by the model.

This capability transforms transformers from mere predictors into scientific inquiry tools. Researchers can analyze attention scores to uncover complex data patterns.

For instance, in protein evolution, a transformer trained to distinguish ancestral and descendant sequences can use attention to identify key amino acid positions. These high-attention residues, deemed biologically salient, can guide lab experiments, highlighting sites critical for functional shifts or adaptation.

## A Practical Example in Keras/TensorFlow

The `retrieve_attention_scores.py` script demonstrates this process:

1.  **Model and Layer Extraction**: Loads a pre-trained transformer and accesses its embedding and encoder layers. This is vital for inspecting internal model workings.
2.  **Manually Forward-Propagating Data**: Input data is passed step-by-step through layers. Crucially, the final encoder is called with `return_attention=True` to output its calculated attention weights.
3.  **Processing the Attention Scores**: Raw attention scores have a shape of `(batch, heads, seq_len, seq_len)`.
    *   Multiple "attention heads" (for different relationship types) are averaged (`np.mean(..., axis=0)`) for a consolidated importance view.
    *   The resulting matrix shows how much attention position `j` pays to position `i`. Averaging over columns (`np.mean(..., axis=0)`) gives a single attention score per residue, representing total attention received.
4.  **Saving and Analysis**: Final per-residue scores are saved to CSV. This allows visualization (e.g., plotting scores against residue position). Peaks indicate amino acids the model found most salient, generating testable hypotheses.

In summary, extracting and analyzing transformer attention scores turns a predictive model into an investigative tool, offering a data-driven method for insights and guiding research.
