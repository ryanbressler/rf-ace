**The manual pages have been written on the basis of RF-ACE version 1.1.0**

# Installation #

Download the latest stable release from the [download page](http://code.google.com/p/rf-ace/downloads/list), or checkout the latest development version (to directory rf-ace/) by typing
```
svn checkout http://rf-ace.googlecode.com/svn/trunk/ rf-ace
```

Compiler makefiles for Linux (`Makefile`) and Visual Studio for Windows (`make.bat`) are provided in the package. In Linux, you can compile the program by typing
```
make
```

In Windows and using Visual Studio, first open up the Visual Studio terminal and execute either `make_win32.bat` or `make_win64.bat` by typing
```
make_win32
```
or
```
make_win64
```
Simple as that!

If you feel lucky, check for compiled binaries at the [download page](http://code.google.com/p/rf-ace/downloads/list).

# Supported data formats #
RF-ACE currently supports two file formats, Annotated Feature Matrix (AFM) and Attribute-Relation File Format (ARFF).

**Annotated Feature Matrix (AFM)**

Annotated Feature Matrix represents the data as a tab-delimited table, where both columns and rows contain headers describing the samples and features. Based on the headers, the AFM reader is able to discern the right orientation (features as rows or columns in the matrix) of the matrix. Namely AFM feature headers must encode whether the feature is (`N`)umerical, (`C`)ategorical, or (`B`)inary, followed by colon and the actual name of the feature as follows:

  * `B:is_alive`
  * `N:age`
  * `C:tumor_grage`

In fact any string, even including colons, spaces, and other special characters, encodes a valid feature name as long as it starts with the preamble `N:`/`C:`/`B:`. Thus, the following is a valid feature header:

  * `N:GEXP:TP53:chr17:123:456`

Sample headers are not constrained, except that they must not contain preambles `N:`/`C:`/`B:`, being reserved for the feature headers.

Missing values can be any upper/lowercase combination of

  * NA
  * NAN
  * NULL
  * ?

And binary/categorical features may have string literal entries. Here is an example of a valid AFM:
```
	S1	S2	S3	S4	S5	S6	S7	S8	S9	S10
N:F1	NA	8.5	3.4	7.2	5	6	7	11	9	NA
N:F2	2	3	4	5	6	NA	NA	9	NA	10
B:F3	NA	NA	NA	NA	MALE	FEMALE	MALE	MALE	MALE	FEMALE
N:F4	10	9.9	8	7	6	5	4	3	2.4	1
C:F5	3	3	3	4	4	5	3	2	2	2
N:F6	9	8	7	9	8	7	3	2	1.0	99.23
```

**Attribute-Relation File Format (ARFF)**

When the feature matrix is as ARFF, the aforementioned constraints for feature headers do no longer apply, since the ARFF itself has an internal format for specifying feature type.

Read more about [ARFF specification](http://www.cs.waikato.ac.nz/~ml/weka/arff.html).

# Parameter print-out #

```
Forest Options:
 -f / --forestType       Forest type: RF (default), GBT, or CART
 -n / --nTrees           [RF and GBT only] Number of trees in the forest
 -m / --mTry             [RF only] Fraction of randomly drawn features per node split
 -a / --nMaxLeaves       Maximum number of leaves per tree
 -s / --nodeSize         Smallest number of train samples per leaf node
 -k / --shrinkage        [GBT only] Shrinkage applied to evolving the residual
 -c / --contrastFraction [Filter only] the fraction of contrast features sampled to approximate the null distribution
File Options:
 -F / --filterData       Load data file (.afm or .arff) for feature selection
 -I / --trainData        Load data file (.afm or .arff) for training a model
 -w / --featureWeights   Load feature weights from file
 -W / --whiteList        Load white list from file
 -B / --blackList        Load black list from file
 -T / --testData         Load data file (.afm or .arff) for testing a model
 -L / --loadForest       Load model from file (.sf)
 -V / --saveForest       Save model to file (.sf)
 -A / --associations     Save associations to file
 -P / --predictions      Save predictions to file
 -R / --pairInteractions Save pair interactions to file
 -G / --log              Save log to file
Filter Options:
 -p / --nPerms           Number of permutations in statistical test
 -t / --pValueTh         P-value threshold in statistical test
 -o / --importanceTh     Importance threshold
General Options:
 -i / --target           Name or index for the variable to make modeling for
 -D / --dataDelim        [AFM only] Field delimiter
 -H / --headDelim        [AFM only] Feature type and name delimiter
 -X / --pruneFeatures    Prune Features
 -S / --seed             Seed for random number generator
 -e / --nThreads         Number of threads if using multithreading
 -R / --maxThreads       Flag to make use of all available threads
 -d / --defaultWeight    Default feature weight, if using feature weighting
```

# Identify important features from data using RF-ACE #
RF-ACE performs feature selection for a given target and feature matrix using `rf-ace` with nTrees=10 and mTry=10 as follows:
```
bin/rf-ace -F data.arff -i target -n 100 -m 10 -A associations.tsv
```
The format of the association list is as follows (one entry per row):
```
TARGET    PREDICTOR    P-VALUE    IMPORTANCE    CORRELATION    NSAMPLES
```
where `NSAMPLES` is the number of non-missing samples shared by the target and predictor. There will be one such entry (row) per each identified, statistically significant association, sorted by p-value to increasing order.

## Giving target feature as an integer ##
It is possible to replace the name of the feature as the argument for `-i / --target` with a `0`-base integer pointing to the proper row/column of the data file. So valid arguments for `-i` are `0,1,...,N-1`, where `N-1` is the number of features in the input data file. The following command would take the first feature in the file as the target feature and regard others as candidate predictors:
```
bin/rf-ace -F data.arff -i 0 -n 100 -m 10 -A associations.tsv
```

## Tuning parameters of the statistical test ##
By default `rf-ace` assigns the p-value threshold to `0.05`, but which can be manually tuned by the `-t / --pValueTh` parameter like so:
```
bin/rf-ace -F data.arff -i target -n 100 -m 10 -A associations.tsv -t 0.001
```
If you would like to have more permutations (default `20`) sampled for the t-test, you can do that using the parameter `-p / --nPerms` like so:
```
bin/rf-ace -F data.arff -i target -n 100 -m 10 -A associations.tsv -p 50
```
In many circumstances, however, it is not recommended to reduce the number of permutations; consider `20` to be the lower limit at which the statistical test behaves robustly across wide range of situations. If RF-ACE recognizes the number of permutations is not sufficient for the t-test, an error will be raised and the program exits; increasing the number of permutations for the t-test will most likely remedy the issue.

## Removing features with too few samples overlapping with the target ##
As the input data can be sparse, and depending on which target has been selected, it may become favorable to remove features that are too sparse for reliable analysis. By default `rf-ace` removes features with less than `5` non-missing samples shared by the feature and target, but which can be set to any other value using the parameter `-X / --pruneFeatures` like so:
```
bin/rf-ace -F data.arff -i target -n 100 -m 10 -A associations.tsv -X 60
```
which would thus cause features with less than `60` intersecting non-missing samples to be removed from analysis.
### Blacklisting features ###
In some cases the input data contains features that for some reason are not suited for as predictors, but whose explicit removal from the input data may be unfavorable. If that is the case, then one can invoke the blacklisting parameter `-B / --blackList` like so:
```
bin/rf-ace -F data.arff -i target -n 100 -m 10 -A associations.tsv -B blacklist.txt
```

# Train a Random Forest (RF) model & predict #
Where the use of the file input option `-F / --filterData` triggers variable selection, `-I / --trainData` triggers model training. In conjuntion with some test data, which can be given using the test input file option `-T / --testData`, predictions can be written to a file using the prediction output option `-P / --predictions` as follows:
```
bin/rf-ace -I data.arff -i target -n 100 -m 10 -T testdata.arff -P predictions.tsv
```
It is also possible to save the predictor for a later use with the save forest option `-V / --saveForest`:
```
bin/rf-ace -I data.arff -i target -n 100 -m 10 -V forest.sf
```
Now that the predictor is saved in a file, which is by the way both human- and machine-readable, it can be loaded for future prediction tasks with the load forest file option `-L / --loadForest` as follows:
```
bin/rf-ace -L forest.sf -T testdata.arff -P predictions.tsv
```

# Train a Gradient Boosting Tree (GBT) model & predict #
RF-ACE uses Random Forest as the default model for training and prediction; the model used can be specified explicitly using the forest type option `-f / --forestType`:
```
bin/rf-ace -I data.arff -i target -f RF -n 100 -m 10 -V forest_RF.sf
```
which trains a Random Forest model, or
```
bin/rf-ace -I data.arff -i target -f GBT -n 100 -a 10 -k 0.01 -V forest_GBT.sf
```
which trains a Gradient Boosting Tree model with `-a / --nodeSize` of 10, and `-k / --shrinkage` of 0.01. Using the trained and saved GBT predictor happens the same way as with RFs:
```
bin/rf-ace -L forest_GBT.sf -T testdata.arff -P predictions.tsv
```