## Overview ##

RF-ACE is a Random Forest [[1](http://oz.berkeley.edu/~breiman/randomforest2001.pdf)] implementation that handles numerical, categorical, and textual variables with missing values.

A main feature of the implementation is variable filtering. Artificial contrasts are used to assess statistical significance of variables. This is where the name ACE comes from, Artificial Contrasts with Ensembles. RF-ACE implements the feature filtering first described in [[2](http://enpub.fulton.asu.edu/workshop/FSDM05-Proceedings.pdf#page=74)], and later in
[[3](http://jmlr.org/papers/v10/tuv09a.html)]. The best variable subset selection described in
[[3](http://jmlr.org/papers/v10/tuv09a.html)] is not implemented.

For regression or classification RF-ACE provides a full conditional distribution of the target variable. This is based on ideas behind Quantile Regression Forests [[4](http://jmlr.org/papers/v7/meinshausen06a.html)], in which the leaf nodes retain the training data samples propagated into them, rather than retaining just the mean or mode of samples. As a result, any quantile of the prediction can be queried.

Text data is preprocessed using the “hashing trick” described in [[5](http://alex.smola.org/papers/2009/Weinbergeretal09.pdf)], except that we use murmurhash3. A one-step stochastic greedy [[6](http://www.thinkmind.org/index.php?view=article&articleid=soft_v4_n12_2011_1)] splitter strategy is used with categorical variables.

Forests can be learned with multiple threads if needed, and saved for later use.

The package also compiles to an R library "rfacer".

[[1](http://oz.berkeley.edu/~breiman/randomforest2001.pdf)] L. Breiman, “Random Forests”, Machine Learning, 45(1):5–32, 2001

[[2](http://enpub.fulton.asu.edu/workshop/FSDM05-Proceedings.pdf#page=74)] Eugene Tuv, Kari Torkkola, “Feature filtering with ensembles using artificial contrasts”, Proceedings of the SIAM 2005 Int. Workshop on Feature Selection for Data Mining, pages 69-71, April 23, 2005

[[3](http://jmlr.org/papers/v10/tuv09a.html)] Eugene Tuv, Alexander Borisov, George Runger, Kari Torkkola, “Feature Selection with Ensembles, Artificial Variables, and Redundancy Elimination”, Journal of Machine Learning Research 10, 1341--1366, 2009.

[[4](http://jmlr.org/papers/v7/meinshausen06a.html)] Nicolai Meinshausen, “Quantile Regression Forests”, Journal of Machine Learning Research 7, 983-999, 2006

[[5](http://alex.smola.org/papers/2009/Weinbergeretal09.pdf)] Kilian Weinberger, Anirban Dasgupta, John Langford, Alex Smola, Josh Attenberg,  "Feature Hashing for Large Scale Multitask Learning", Proc. ICML 2009

[[6](http://www.thinkmind.org/index.php?view=article&articleid=soft_v4_n12_2011_1)] Viswanathan, Viswa, Anup Sen, and Soumyakanti Chakraborty. "Stochastic Greedy Algorithms: A leaning based approach to combinatorial optimization." International Journal On Advances in Software 4.1 and 2 (2011): 1-11.

## Installation ##
To get started in Linux, grab the latest source and `make`:
```
svn checkout http://rf-ace.googlecode.com/svn/trunk rfacer
cd rfacer/
make
```
Sometimes compiling with threading fails, in which case you can use
```
make no-threads
```
The source code also compiles into an R package (after installing [Rcpp](http://cran.r-project.org/web/packages/Rcpp/index.html) from CRAN):
```
R CMD INSTALL rfacer
```
For more information see the documentation [here](http://code.google.com/p/rf-ace/wiki/Manual) and [here](https://code.google.com/p/rf-ace/wiki/RPackageInstallation).

### Issues with threads? ###
[Compile and Run Branch "thread-fixes"](http://code.google.com/p/rf-ace/wiki/CompileThreadFixes)
Instructions have been posted to compile the "thread-fixes" branch of version 1.2.5 for Mac OS X and Linux.  The branch does not support Windows or the Rcpp module.

```
svn checkout http://rf-ace.googlecode.com/branches/thread-fixes rfacer
cd rfacer/
# follow instructions on https://code.google.com/p/rf-ace/wiki/CompileThreadFixes
```


### Join the Google Groups [mailing list](http://groups.google.com/group/rf-ace) to receive updates by e-mail ###

## Case study: Large-scale data exploration in The Cancer Genome Atlas ##
In a joint effort together with **[Tampere University of Technology](http://www.tut.fi)** and **[Institute for Systems Biology](http://www.systemsbiology.org/)**, associations uncovered within **[The Cancer Genome Atlas](http://cancergenome.nih.gov/)** using RF-ACE can be viewed at **[Regulome Explorer](http://explorer.cancerregulome.org)**, an interactive web application developed for exploring associations across molecular features spanning the human genome. With help of **[Techila](http://www.techila.fi)** and **[Golem](http://code.google.com/p/golem/)**, CPU intensive but embarassingly parallel computation was distributed across a collection of ~1000 CPUs, cutting down computation from years to days.

## 28th June 2012: RF-ACE featured at Keynote by Urs Hölzle in Google I/O 2012 ##
_During that Google I/O demo, by using 600,000 cores, that genomics app was probably hitting the Compute Engine for somewhere between 1 and 10 petaflops of computational power — between 1 and 10 quadrillion calculations per second. For a few moments, the Compute Engine was probably the third fastest computer in the world._ --[TechCrunch](http://www.extremetech.com/extreme/131962-google-compute-engine-for-2-millionday-your-company-can-run-the-third-fastest-supercomputer-in-the-world)

[More about the case study at GCE developers page](https://developers.google.com/compute/io)

[Press release of the case study](https://cloud.google.com/files/ComputeISBCaseStudy.pdf)




---

## Release history ##
**RF-ACE v1.2.5 release -- Mar 20th 2013:**
  * Quantile Random Forest replaces Random Forest for both classification and regression
  * Predictions take a probabilistic form
    * For classification, class probabilities are estimated
    * For regression, quantiles are estimated
      * Optionally, samples used for generating the quantiles can be output to file
  * Forest file (.sf) now contains information about the decrease of impurity at each split point
  * Gradient Boosting Trees will be temporarily out of use
  * Variable importance now uses mean decrease in impurity (MDI) instead of mean minimal depth (MMD).
  * R package extended to support latest changes

**RF-ACE v1.1.0 release -- Dec 5th 2012:**
  * RF-ACE can be used as a library
    * R wrapper provided as a package (will be submitted as CRAN package)
  * Implemented mean minimal depth score as a robust alternative for variable importance
    * An option to switch between importance and mean minimal depth will be added
  * Redesign of command line arguments ( see help/manual for details )
  * Lots of refactoring

**RF-ACE v1.0.7 release -- Aug 28th 2012:**
  * Fixed one naming convention issue between standard libraries that ship with different GCC versions
  * T-test uses Welch-Satterthwaite's approximation to account for unequal population variances
  * Implemented switch to not launch any child threads if nThreads == 1
    * saves resources
  * NOTE: Currently compiles nicely at least in Linux and with GCC >= 4.4.6 !!
    * Will try to get MinGW working under Windows

**RF-ACE v1.0.6 release -- Aug 17th 2012:**
  * Migrated to C++0x standard. Making use of
    * hash tables in string->int lookups
    * multithreading in forest growing and prediction
    * random number generators => third-party implementations obsolete
    * special functions => third-party implementations obsolete
  * Corrected the way RF-ACE aggregates trees into a single predictor
  * Contrast summary now mean instead of median
  * Removed calculation of variance in split point selection
    * Eliminates overflows with large sample sizes
  * Regular CART (single tree) can now be grown
    * See help for more info
  * RF now grows unpruned trees to (almost) maximal extent
  * Benjamini-Hochberg p-value correction implemented as an option


**RF-ACE v1.0.4 release -- March 25th 2012:**
  * Building either RF ( -R / --RF ) or GBT ( -G / --GBT ) predictor now possible
    * Example: bin/rf-ace-build-predictor --RF ...
    * Default is RF
  * Optimized numerical splitter for numerical targets (~30% boost)
  * Revised background sampling for statistical testing
  * Complete re-implementation of the t-test
    * Previous version was making some false assumptions
  * Takes tied samples into account when testing for splitting
  * Predictor builder prints out OOB error and TOTAL error
    * Good for assessing how well the predictor generalizes to new data
    * In practice
      * With better-than-random predictors: OOB error < TOTAL error
      * With random predictors: OOB error ~ TOTAL error
      * With overfitting predictors: OOB error > TOTAL error
  * Changes in user interface parameters
    * NEW: seed ( -S / --seed ) for the random number generator (Mersenne Twister)
      * By default seed == system clock + elapsed CPU cycles
    * RF and GBT parameters share the same names
    * mTry is specified as positive integer
  * Usage examples are printed when help becomes invoked
  * When all features become pruned, log is updated and program exits normally
  * Lots of re-factoring of code

**RF-ACE v1.0 release -- February 14th 2012:**
  * Optimizations and default parameter tweaks: **over 10x speed-up!**
  * Blacklists/whitelists are now working with index list inputs
  * Three new programs:
    * rf-ace-filter
    * rf-ace-build-predictor
    * rf-ace-predict
  * Updated forest file (.sf) format
    * Generalizes better to various feature naming etc. conventions
  * Create and save predictor (.sf) with “rf-ace-build-predictor”
  * Load forest predictor (.sf) with “rf-ace-predict” and make predictions with novel data
  * Replaced exhaustive search of a binary split with categorical splitter with a greedy one
    * Computational complexity linear as a function of cardinality of the splitter !
  * Better factoring of code with new namespaces
  * Fixed a bug that resulted in incorrect interpretation of mathematical expressions by Visual Studio compiler

**RF-ACE v0.9.9 release -- February 2nd 2012:**
  * Binaries for x86 and x64 Windows XP available
  * Updated forest output writer ( -F / --forest )
    * The format for the forest is almost stable; once fully stable, I will implement a forest reader, meaning one can then save the built model for, say, prediction
  * New whitelist/blacklist functionality with which to fine-tune predictor space
    * -W / --whitelist and -B / --blacklist
  * Updated ARFF reader seg. fault bug
    * This was readily fixed in the intermediate version 0.9.8b
  * Updated log outputs ( -L / --log )
  * Refactored namespaces and classes (namespace math, class statistics::RF\_statistics)
    * This is still far from finished
  * Updated node counter
    * At least works faster, but may still contain bugs (nothing dramatic, though)
  * By default features with less than 5 shared samples with the target will be removed
    * Tune this manually with -X / --prune\_features
  * Fixed a small bug when percolating samples through the trees
    * This affected calculation of importance scores in the presence of missing values
    * Expect to see more association on the outputs!
  * Simplifications in t-test implementation, but no functional change
  * Fixed one implicit type cast, which was causing seg. fault in 64bit Windows
  * Extended manual page to include an explicit example of the AFM data format

**RF-ACE v0.9.8 release -- January 10th 2012:**
  * Possibility to predict novel measurements ( -T / --testdata )
  * Possibility to turn feature selection with RFs off ( -N / --noFilter )
  * Reduced redundancy in class interfaces
  * Refactored main program
  * Eliminated a bug while parsing ARFF data
  * default p-value threshold changed from 0.1 to 0.05
  * GBT forest print-outs ( -F / --forest ) updated to include the name of the target feature
  * The log file ( -L / --log ) now includes:
    * RF-ACE version
    * mean and std importance score for real features
    * mean and std importance score for contrast features
    * mean number of nodes per tree
    * nodes created per second

**RF-ACE v0.9.7 release -- December 29th 2011:**
  * Several small updates, aiming to make prediction with novel data more straightforward, has been made.
  * Data prediction with novel data has been turned OFF until all planned updates are in-place
  * By applying a no-filter flag ( -N / --noFilter ), feature filtering with RFs can be turned OFF (GBT will be used with all available features)
  * By providing a forest output file ( -F / --forest ) the GBT forest will be written to the file
  * By providing a log output file ( -L / --log ) the log will be written to the file; log is currently empty
  * It is now much easier to specify which outputs the user wants (associations, GBT predictor, predictions, log, etc.)

**RF-ACE v0.9.6 release -- December 18th 2011:**
  * The fraction of sampled candidate contrast features for splitting is now tuned from 50% down to just 1%. This will guarantee significantly better trees being grown, while keeping the fraction of contrasts in the trees large enough so that the null distribution can reliably be constructed. You will notice that for significant associations the p-values are now much closer to zero, indicating increased ability to separate true signal from noise.
  * Data and header delimiters can now be changed ( '\t' and ':', respectively, are the defaults )
  * Started implementing a logging feature, with which quality of the analysis can easily be assessed.
  * The default number of candidate features for splitting, mTry, is now 10% of the number of features, i.e. significantly more than it used to be ( sqrt(nFeatures) ). This also guarantees better trees, however:
  * The algorithm now runs slower due to increased CPU load. I will concentrate on cutting down the computing time once I get the planned updates finished.

**RF-ACE v0.9.5 release -- November 14th 2011:**
  * Killed a bug (split decision was 100% biased towards the "left" leaf when the data point was NA; the fix now assigns the sample randomly to left and right according to the fraction of training samples in left and right, respectively) that was severely degrading the prediction accuracy of the algorithm
  * consequently, killing that bug also improved the accuracy of identifying associations

**RF-ACE v0.9.3 release -- November 8th 2011:**
The algorithm has reached a level that I'm fairly confident to say it's nearly bug-free and contains almost all the features planned to be implemented. The latest update features some essential bug fixes:
  * node splitting with numerical features sometimes resulted in under-indexing a vector
  * categorical splitter was calculating the split fitness partly wrong

Both of these bugs, which are now eliminated, were severely decreasing the qualities of the trees they occurred in, but as RF-ACE is a tree ensemble learning algorithm, thus consisting of thousands of trees, the average performance wasn't affected much.

One larger, yet missing feature is support for multiple splitters (main splitter + surrogates) per tree junction, which is supposed to yield better performance with highly sparse data. Support for surrogates will be added in near future.

**RF-ACE v0.9.1 release -- November 6th 2011:**
  * Splitting with features is now implemented as it is in the original formulation
    * if the user wants, (faster) split approximation can be turned ON
  * Gray Code implementation for efficient split testing with categorical features
  * Updated default parameters for the Random Forest

**RF-ACE v0.8.5 release -- October 5th 2011**:
  * Now RF-ACE runs exactly once, for a uniquely specified target
  * Default test has been changed back to the t-test
  * Default number of permutations is 20
  * Improved handling of NaN's and negative importance scores in statistical testing
  * Modified print-outs
  * Lots of small tweaks

**RF-ACE v0.8.0 release -- September 14th 2011**:
  * Restructured the logic of the main program
  * Data prediction works better now
  * Improved print-outs
  * Updated help

**RF-ACE v0.7.5 release -- August 29th 2011**:
  * GBT is now functional in data prediction, so yes:
  * RF-ACE predicts with new data
  * Unit testing is introduced, making development more organized
  * Tons of small updates and bug fixes
  * Started working on support for feature masks (for exclusion of features from analysis)

**RF-ACE v0.5.5 release -- July 5th 2011**: much has changed since the last version:
  * Node class is now dynamic, making tree construction smoother and more memory-efficient
  * GBT is now part of the main program of RF-ACE, albeit not fully functional yet
  * As target one can now specify a string that will be grepped with feature headers
    * if multiple feature headers match the string, multiple RF-ACE calls are made, and results concatenated to the specified output file
  * Fixed a bug that was annoyingly making contrasts to never enter the trees
  * Lots of small tweaks

**RF-ACE v0.4.0 release -- July 1st 2011**: The next stable release of RF-ACE, version 0.4.0, is ready. Although functionally very similar to v0.3.5, most of the internal components have been revised and simplified, and naming conventions unified. There will be a few more of such revisions, this time concentrating on making tree generation more dynamic and shifting splitting functions under the Node implementation.

**June 29th 2011**: a major revision of the internal structure RF-ACE is now completed, and code has been committed to the trunk. Also, a makefile for Visual Studio command line compiler (cl) is provided. Some further updates will be executed before the next stable release will be announced. Stay tuned.

**RF-ACE v0.3.5 release -- June 24th 2011**: version 0.3.5 is out! Lots of small tweaks since the last version:
  * fixed a bug that allowed the target feature to enter the list of predictors
    * the target itself didn't exist anywhere in the trees
  * sufficient nodesize is now estimated directly from the data
    * the algorithm now adapts to larger sample sizes by tuning nodesize up
    * as nodes grow bigger, trees become smaller and leave room for more permutations
  * increased default permutation size from 20 to 50
    * increases statistical power
  * finished implementing ARFF support, it should be working now
  * reformatted print-outs


**RF-ACE v0.3.0 release -- June 21st 2011**: RF-ACE has now reached version 0.3.0 (check the source package and Win32 binary [here](http://code.google.com/p/rf-ace/downloads/list)). It can be considered a stable release of the algorithm. RF-ACE v0.3.0 does the following:
  * accepts AFM (Annotated Feature Matrix) files as inputs
    * identifies the type of the input automatically
  * identifies the orientation of the input, if AFM, automatically
  * handles numerical and categorical features
  * handles missing values
  * estimates ntrees and mtry based on data dimensions and the number of missing values
  * constructs multiple Random Forests
    * an optimized implementation of the original RF (less sorting involved)
  * uncovers statistically significant associations using Mann-Whitney U-test