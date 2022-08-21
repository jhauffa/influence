# Social Influence Experiments

This repository contains the source code of the main experiments conducted in the course of a PhD thesis:

(A link to the thesis will be inserted here, pending its official approval.)

The first two experiments approach the problem of social influence detection at different levels of granularity: At the micro level, we attempt to test if a change in the behavior of one person can be explained by the earlier behavior of another, which would be evidence for a causal relationship in the sense of Granger causality. At the meso level, the detection of social influence can be reframed as a prediction problem: Which parts of a person's social environment are most useful for predicting that person's future behavior? The third and final experiment deals with the temporal distances of interactions on social media.

## Contributors

The following people, listed in chronological order, were involved in the development:

* Jan Hauffa (jhauffa@gmail.com)
    - research and development, maintenance of the codebase
* Benjamin Koster
    - topic modeling, meso-level social influence experiment
* Florian Hartl
    - tokenization, topic modeling, meso-level social influence experiment
* Valeria Kalteis
    - meso-level social influence experiment
* Julia Strauß
    - maximal cliques/plexes, meso-level social influence experiment
* Felix Sonntag
    - meso-level social influence experiment

## Publications

* Jan Hauffa, Benjamin Koster, Florian Hartl, Valeria Köllhofer, Georg Groh: "Mining Twitter for an Explanatory Model of Social Influence", 2016, 2nd International Workshop on Social Influence Analysis @ IJCAI
* Jan Hauffa, Wolfgang Bräu, Georg Groh: "Detection of Topical Influence in Social Networks via Granger-Causal Inference: A Twitter Case Study", 2019, Workshop on Social Influence @ ASONAM
* Jan Hauffa, Georg Groh: "A Comparative Temporal Analysis of User-Content-Interaction in Social Media", 2019, 5th International Workshop on Social Media World Sensors @ HT

The meso-level social influence experiment is described in the 2016 paper, and the micro-level influence experiment in the 2019 ASONAM paper.

## Usage

This codebase is mainly intended to be used as a reference, to be read alongside the PhD thesis, which cannot possibly describe all implementation details. If you want to reproduce the results locally, or run the experiments on a dataset of your own, follow the instructions below.

### Prerequisites

1. Run `mvn package` to build a JAR file `target/influence-1.0-SNAPSHOT.jar` that contains the compiled code and all its dependencies.
2. Provide the data you wish to analyze in a MySQL database and create a MySQL user account that has read and write access to that database. The [crawler](https://github.com/jhauffa/crawler) repository contains a number of Java programs that generate databases in a suitable format. If you prefer to write your own exporter, you can use class `edu.tum.cs.postprocessing.SocialMediaDao` directly or as a template.
3. Create a file `experiment.properties` that holds the configurable experiment parameters (in Java "key=value" format). Property "DatabaseAccessor.URL" should be set to the JDBC URL of the database. If you are working with Twitter data, set "SocialMediaDaoFactory.dataSource" to "twitter", "TwitterDao.numUsers" to the number of users you want to analyze, and "TwitterDao.usersList" to the name of a text file that holds the user IDs (one ID per line). For any other data source, set "SocialMediaDaoFactory.dataSource" to the common prefix of the names of the tables generated in the previous step.

You can now run any program inside the JAR as follows (on a single line):
`java -Dedu.tum.cs.util.config="path/to/experiment.properties" -cp target/influence-1.0-SNAPSHOT.jar edu.tum.cs.xyz.ClassName`

If you intend to run the meso/micro level influence experiments:

1. Configure the following additional parameters in `experiment.properties`: Set "endDate" to one day past the last day of the observation period (yyyy.MM.dd, e.g., "2012.06.01"). Parameter "IntervalCalculator.numIntervals" specifies the number of intervals into which the observation period is subdivided, and "IntervalCalculator.weeksPerInterval" the length of each interval. Parameter "timePeriods" should be set to a comma-separated list of interval lengths; each value specifies how many days of an interval are actually used. Set "TwitterTokenizer.dataPath" to the directory that contains the list of dataset-specific stop words (should be called `stopWords-en.txt` and contain one word per line). You can use `edu.tum.cs.nlp.tokenizer.GenerateStopWordList` to generate a list of candidate words. For replication of the experiments in the publications listed above, set "tokenizeURLs", "TwitterMessageLoader.resolveURLs", and "TwitterMessageLoader.augmentTweets", and "TopicModel.online" to false. Set "TopicModel.modelPath" to the path of an existing directory, and "TopicModel.numTopics" to a reasonable number of topics (e.g., 150). Note that some parts of the code assume that the number of topics does not exceed 255. If you are working with non-Twitter data, set parameter "GenericMessageLoader.minRecipientGroupSize" to the minimum number of recipients for a message to be considered non-addressive (e.g., 10).
2. Run `edu.tum.cs.nlp.topic.FitART` to build an ART topic model. Depending on the size of the dataset, you may need to increase the JVM heap size or reduce the number of threads via parameter "numThreads".

### Meso-Level Influence Experiment

1. Run `edu.tum.cs.graph.GenerateSocialEdgesInDB` to extract the explicit (if present) and implicit social network graphs. If parameter "GenerateSocialEdgesInDB.alwaysUseCommGraphMetrics" is false, graph metrics are derived from the explicit graph if available. Set to false or leave unset when replicating the 2016 paper.
2. To find maximal cliques in the social network graphs, run `edu.tum.cs.graph.clique.FindSubgraphs`. For replicating the 2016 paper, pass "bidir" as a command line argument (converts a directed graph to an undirected graph by inserting an undirected edge if and only if two nodes are connected bidirectionally). Pass "comm" to analyze the implicit graph. Afterwards, run `edu.tum.cs.graph.clique.FindPercolatedCommunities` to apply the clique percolation algorithm to the set of maximal cliques. Pass the name of the output file generated by `FindSubgraphs` as a command line argument. All maximal cliques have to be loaded into RAM, so you may need to increase the heap size. Finally, run `edu.tum.cs.graph.clique.FindEdgeCommunities` to find communities via edge clustering. Command line arguments "comm" and "bidir" have the same function as in the case of `FindSubgraphs`.
3. If you have sufficient data points (i.e., nodes and edges of the social network graph) for ANOVA, which requires partitioning the sample into non-overlapping subsets for each combination of factor levels, set "InfluenceExperimentsAnova.splitAnova" to true. Then, run `edu.tum.cs.influence.meso.InfluenceExperimentsAnova` to perform the actual influence experiment.
4. Run `edu.tum.cs.influence.meso.eval.AggregateResults` to generate a summary of the results (open in LibreOffice Calc and use "AutoFilter"!) and to identify combinations of experiment parameters that perform worse than the baselines.
5. Run `edu.tum.cs.influence.meso.eval.PrepareAnova` to convert the results to a suitable format for ANOVA. Pass `discarded-variants.csv` (generated by `AggregateResults`) as a command line parameter. Unpack the JAR file to obtain the R script `visualization/anova.R`. If you set "InfluenceExperimentsAnova.splitAnova" to true in step 3, you can run the script directly. Otherwise, comment out the lines starting with `processExperimentResults(` at the end of the file, uncomment the lines that start with `processExperimentResultsSimple(`, and adjust the parameters to match your data.

### Micro-Level Influence Experiment

1. Set parameter "InfluenceMeasurement.mutualNetworksOnly" to true to only analyze the undirected social network graphs formed by reciprocal communication / explicit edges. Set "GrangerInfluenceModel.absMinObservations" to the minimum number of observations that have to be present for an edge to be considered for analysis (e.g., 6). Similarly, "GrangerInfluenceModel.relMinObservations" can be set to a value > 0 to specify a  minimum proportion of present edges. Set "InfluenceMeasurement.noiseTopicThreshold" to a value > 0 to ignore topics if their corresponding components of prior parameter vector alpha are below that threshold.
2. Run `edu.tum.cs.influence.micro.InfluenceMeasurement` to perform the experiment.
3. Run `edu.tum.cs.influence.micro.ExtractTweetSample` to extract the tweets of random user pairs within the observation period and to generate an annotation sheet into which the results of the manual assessment can be entered. After annotation is completed, create separate files for rows with positive and negative annotations.
4. Run `edu.tum.cs.influence.micro.EvaluateNetworks` with the negative annotation file as the first and the positive annotation file as the second command line argument. If there are no positive annotations, the second argument can be omitted.

### Temporal Analysis of User-Content-Interaction

1. To extract sequences of temporal distances between interactions from a social media dataset, run `edu.tum.cs.time.hmm.BuildInteractionChains`. The first command line argument is the prefix of the output files, the second argument identifies the data source (same as parameter "SocialMediaDaoFactory.dataSource"), the third argument is the interaction type (one of REPLY, SHARE, BOTH).
2. Run `edu.tum.cs.time.hmm.FitInteractionHmm` to fit an HMM to the sequences. First command line argument is the string "fit", second argument is the name of the sequence file generated in the previous step.
3. Finally, run `edu.tum.cs.time.hmm.HmmPlotter` to plot and analyze the HMM. The program takes the name of the HMM parameter file generated in the previous step (typically "hmm-param.bin") as an argument.
