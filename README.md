## Methods

We applied a machine learning approach to identify genes associated to response to treatment, and to find a list of genes that could be potential predictors of response.  
### Data preprocessing. 
The microarray data was prepared with Asa pipeline, bulk RNAseq data was log2 transformed (+1 to avoid infinites in the log2).
Both microarray dataset and bulk dataset samples were median centered and sd scaled. 
To reduce the risk of overfitting the model, we have crossvalidate each of the following step.
### Differential gene expression analyis. 
As an first exploratory analysis, and as a feature selection before the machine learning models, we have run a differential gene expression analysis.   
We have run a GLM model with the limma package inside a cross validation loop. The data was split in split in 10 folds (caret::createFolds, the random sampling is done within the levels of y when y is a factor in an attempt to balance the class distributions within the splits).  
The we run all the combination of 8 out of 10 (45 runs). Finally, we compute intersect (max pvalue), union (min pvalue), and average pvalue of all the runs.   
We included the batch as a coveriate in the limma formula ~ batch + response  
As an high confidence DE results we have the intersection (gene significative in all runs).  
As a less stringent selection, for feature to input to the machine learning models, we take average pvalue (not corrected) threashold: 0.05   

