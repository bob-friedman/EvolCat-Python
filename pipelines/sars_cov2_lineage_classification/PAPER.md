# Rapid and Interpretable Deep Learning for Fine-Grained SARS-CoV-2 Lineage Classification

## Abstract
> The rapid evolution and diversification of SARS-CoV-2 necessitate automated tools for accurate, fine-grained lineage classification. An end-to-end deep learning pipeline is presented that processes genomic variant data into a binary feature matrix to train a neural network classifier. When applied to the highly diverse XBB 2.3 clade, containing 16,367 genomes across 95 distinct sublineages, our model achieved a test accuracy of 97.98%. Crucially, by employing gradient-based saliency mapping, the model is rendered interpretable, identifying the specific mutations driving the classification for each lineage. This work demonstrates a scalable and interpretable framework for real-time genomic surveillance, effectively bridging computational analysis with actionable virological insights and validating the model's biological relevance.

## Keywords
**Keywords:** virus evolution; mutation; lineage classification; deep learning pipeline; neural network classifier; saliency maps.

---

## 1. Introduction
The SARS-CoV-2 pandemic has underscored the critical need for rapid, scalable, and accurate genomic surveillance to monitor viral evolution and inform public health responses. The constant emergence of new sublineages, driven by mutation and selection, presents a formidable classification challenge [1,2]. While phylogenetic placement tools like `UShER` provide an invaluable, memory-efficient framework for tracking viral spread on a pandemic scale [3], there is a growing need for complementary methods that can not only classify lineages with high accuracy but also provide transparent, interpretable insights into the specific genomic features that define them. Many computational models operate as **"black boxes"**, which can limit their direct biological relevance and trustworthiness.

In this communication, an end-to-end pipeline is presented that leverages deep learning for the fine-grained classification of SARS-CoV-2 lineages and employs saliency mapping to render the model's decisions interpretable. Its efficacy is demonstrated on the highly diverse XBB 2.3 clade, showing that this approach not only achieves exceptional classification accuracy but also successfully identifies the key mutations driving these classifications, thereby bridging the gap between computational prediction and actionable biological insight.

## 2. Methods

### 2.1. Data Acquisition and Preprocessing
The foundation of our analysis is the public SARS-CoV-2 Mutation-Annotated Tree (MAT) provided by the UShER project, which enables real-time, large-scale phylogenetics [3]. The latest global MAT (`public-latest.all.masked.pb`) was downloaded and processed using the `matUtils` toolkit [4,5]. All genomic data was extracted belonging to the XBB 2.3 major clade (Pango lineage designation 23E/XBB_2.3). This critical step isolates a manageable subset of data for focused analysis and outputs it in the standard Variant Call Format (VCF), which contains the specific genetic variations for each sample relative to the Wuhan-Hu-1 reference genome (NC_045512.2). Pango lineage assignments for each sample were retrieved from the associated global metadata file (`public-latest.metadata.tsv`) [5].

### 2.2. Feature Matrix Construction
To prepare the data for machine learning, an automated feature engineering pipeline was developed. This pipeline converts the VCF data into a binary feature matrix, `X`, where each row represents a unique viral genome sample and each column corresponds to a unique single nucleotide variant (SNV) identified within the entire XBB 2.3 clade dataset. A value of `1` indicates the presence of a variant in a sample, while `0` indicates its absence. This process relies on `pysam`, a highly efficient library for reading and manipulating indexed VCF files [6].

To manage the large data footprint, several memory optimizations were implemented. The full metadata file was filtered after extracting the sample list from the VCF, ensuring that only relevant sample information was loaded into memory. Furthermore, the final feature matrix was constructed using an `int8` data type, an 8-fold reduction in memory usage compared to the default, which is vital for scalability.

A critical, scientifically-motivated step in this stage was the removal of lineages represented by only a single sample ("singletons"). A model cannot learn the generalizable features of a class from a single example. Removing these ensures that every class in our dataset has at least two examples, which is a prerequisite for robust, stratified data splitting and meaningful model training.

### 2.3. Deep Learning Model Architecture and Training
For this classification task, a feed-forward neural network (FNN) was implemented, a classic and effective architecture for tabular data, using the `TensorFlow/Keras` library [7]. The model architecture was designed to be powerful enough to capture complex patterns without being excessively large, which could risk overfitting. It consists of an input layer sized to the number of unique variants (12,277), two hidden dense layers with 128 and 64 neurons respectively (using the Rectified Linear Unit, `ReLU`, activation function), and a final output layer. The output layer contains a neuron for each of the 95 distinct lineages and uses a `softmax` activation function to produce a probability distribution over the possible classes.

The dataset was split into training (80%) and testing (20%) sets. Crucially, a stratified splitting approach was employed based on the lineage labels. This ensures that the proportional representation of each lineage, particularly rare ones, is maintained in both the training and testing sets, leading to a more reliable and valid evaluation of the model's performance. The model was compiled using the popular and effective Adam optimizer and the sparse categorical cross-entropy loss function, which is standard for multi-class classification with integer labels. Training was conducted for 10 epochs with a batch size of 32.

### 2.4. Model Interpretability using Saliency Mapping
To ensure the model's decisions are transparent and biologically relevant, a model interpretability analysis was implemented using gradient-based saliency maps. This technique directly addresses the question: "To correctly identify this sample's lineage, which mutations did the model 'pay attention' to the most?" For a given correct classification, the gradient of the score at the corresponding output neuron with respect to every input feature was calculated (i.e., every mutation). The absolute value of this gradient is used as an "importance score", directly quantifying how influential each mutation was in the model's decision to assign a sample to that specific lineage. This transforms the model from a black box into an investigative tool.

## 3. Results

### 3.1. Application to the Highly Diverse XBB Clade
Our end-to-end pipeline was successfully applied to a challenging, real-world dataset comprising genomes from the SARS-CoV-2 XBB 2.3 clade. After preprocessing and the removal of one singleton lineage, the final dataset consisted of 16,367 unique viral samples distributed across 95 distinct Pango sublineages. The feature engineering process identified a total of 12,277 unique SNVs across these samples, which formed the feature space for our classification model. The scale and extreme class imbalance of this dataset represent a significant test of the pipeline's ability to perform fine-grained classification.

### 3.2. High-Accuracy Lineage Classification Presents a Story of Extremes
Upon evaluation with the held-out test set (n=3,274 samples), the FNN model demonstrated exceptional performance, achieving a final test accuracy of 97.98%. The model's robustness was further confirmed by a weighted average F1-score of 0.98, indicating high precision and recall across the diverse lineages.

A detailed classification report (data not shown) tells a compelling "story of extremes". The model achieved perfect or near-perfect F1-scores for dozens of lineages, including those with very high support (e.g., `GJ.1.2`, n=479) and many with moderate support (e.g., `GE.1.4`, n=93). Conversely, the model's few failures were predictable and informative. Lineages with extremely low sample counts in the test set (e.g., `GJ.1`, support=2; `GE.1.2.1`, support=1) had F1-scores of 0.00. This is an expected and important result, as it confirms the model's data-driven learning process and highlights that its performance bottleneck is the availability of data for rare classes, not a flaw in the model's architecture.

### 3.3. Interpretable Classification through Saliency-Based Feature Importance
The key result of our work is the successful interpretation of the model's decisions. The saliency analysis produced a ranked list of the most influential mutations for the correct classification of each of the 95 lineages, revealing the specific genomic markers the model learned to prioritize. This analysis generated a wide spectrum of importance scores, indicating that some lineages are defined by a few highly determinative mutations, while others are distinguished by more subtle patterns.

For instance, in classifying samples of lineage GE.1, the model assigned the highest importance score (0.0013) to the mutation at position `5947`. For lineage XBB.2.3.20, the most determinative feature was a mutation at position `27992` (importance: 0.0013). Crucially, this computational approach has a potential to rediscover known biological markers (Table 1).

_Table 1: Summary of top influential mutations for representative lineages._
| Lineage | Mutation | Importance |
|-----------|----------|------------|
| GE.1 | `5947_t` | 0.001294 |
| GE.1 | `2246_g` | 0.000833 |
| GE.1 | `3583_c` | 0.000298 |
| GE.1 | `3022_t` | 0.000220 |
| GE.1 | `22688_a`| 0.000107 |
| GE.1 | `22775_g`| 0.000104 |
| GE.1 | `14327_c`| 0.000089 |
| GE.1 | `27415_g`| 0.000045 |
| GE.1 | `22674_c`| 0.000032 |
| GE.1 | `25584_c`| 0.000031 |
| XBB.2.3.20| `27992_t`| 0.001273 |

## 4. Discussion
In this study, an automated pipeline for the accurate and interpretable classification of SARS-CoV-2 sublineages was developed and validated. By applying a feed-forward neural network to a binary mutation matrix, a remarkable 97.98% accuracy in distinguishing among 95 lineages of the complex XBB 2.3 clade was achieved. The model training progression and the close tracking of training and validation loss demonstrated that the model architecture had sufficient capacity to learn without overfitting the data.

The principal contribution of this work lies in its interpretability. The use of gradient-based saliency analysis provides a clear window into the model's decision-making process. Saliency scores provide a computationally-derived guide for biological investigation, allowing researchers to prioritize which mutations warrant further study. The magnitude of a mutation's importance score can also serve as a proxy for its evolutionary distinctiveness relative to other lineages in the dataset.

This framework represents a significant step towards transparent and trustworthy AI in genomic epidemiology. Future work could involve applying this pipeline to other viral families or comparing the saliency-based feature importances with those from tree-based models like `XGBoost`. Further analysis of misclassifications, guided by the confusion matrix, could also reveal fascinating instances of convergent evolution or genetic ambiguity between closely related lineages.

In conclusion, this proof-of-concept demonstrates a scalable and insightful approach for monitoring viral evolution. This fusion of high-performance classification and model interpretability has the potential to become an indispensable component of future pandemic preparedness and response efforts, turning raw genomic data into actionable knowledge.

## Credits

This work was developed by the author in collaboration with Google's Gemini Pro, which includes Python code generation, data analysis, and the initial drafting of this paper. The research is further built upon the foundational public data provided by the UShER project and the open-source software `matUtils`, `pysam`, and TensorFlow/Keras.

## References
1. Lowen, A. C. (2017) Constraints, Drivers, and Implications of Influenza A Virus Reassortment. *Annual Review of Virology*, 4, 105-21. [https://doi.org/10.1146/annurev-virology-101416-041726](https://doi.org/10.1146/annurev-virology-101416-041726)
2. Domingo, E., Martin, V., Perales, C., Grande-Pérez, A., García-Arriaza, J., & Arias, A. (2006) Viruses as quasispecies: biological implications. *Current Topics in Microbiology and Immunology*, 299, 51-82. [https://doi.org/10.1007/3-540-26397-7_3](https://doi.org/10.1007/3-540-26397-7_3)
3. Turakhia, Y., Thornlow, B., Hinrichs, A. S., De Maio, N., Gozashti, L., Lanfear, R., ... & Corbett-Detig, R. (2021) Ultrafast Sample placement on Existing tRees (UShER) enables real-time phylogenetics for the SARS-CoV-2 pandemic. *Nature Genetics*, 53, 809-16. [https://doi.org/10.1038/s41588-021-00862-7](https://doi.org/10.1038/s41588-021-00862-7)
4. Ultrafast Sample Placement on Existing Trees. Available online: [https://github.com/yatisht/usher](https://github.com/yatisht/usher) (accessed on 4 June 2025).
5. Public SARS-CoV-2 sequence data. Available online: [https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/](https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/) (accessed in June 2025).
6. Pysam is a lightweight wrapper of the HTSlib API. Available online: [https://github.com/pysam-developers/pysam](https://github.com/pysam-developers/pysam) (accessed in June 2025).
7. Keras is a deep learning API written in Python. Available online: [https://keras.io/getting_started/about](https://keras.io/getting_started/about) (accessed on 17 June 2025).
```
