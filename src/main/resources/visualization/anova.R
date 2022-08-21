options(width = 10000)
options(max.print = 1000000)

library(car)
library(emmeans)
suppressMessages(library(ggplot2))

readSCIMData <- function(dir, inputFile, ignoreBaseline) {
	SCIM <- read.table(paste(dir, inputFile, sep = ""),
		header = TRUE, sep = ";", dec = ".",
		colClasses = c("factor", "integer", "factor", "factor", "factor", "factor",
			"factor", "factor", "numeric", "NULL", "NULL", "NULL", "NULL"))
	SCIM$timePeriodBase <- as.ordered(SCIM$timePeriodBase)
	SCIM$timePeriodLength <- as.ordered(as.factor(SCIM$timePeriodLength))

	if (ignoreBaseline) {
		SCIM <- SCIM[which(!startsWith(as.character(SCIM$type), "BASELINE")),]
	}

	# merge type,indicator,weight into single factor; same for fromUserId,toUserId
	SCIM <- within(SCIM, indicatorWeight <- paste(type, indicator, weight, sep = "_"))
	SCIM <- within(SCIM, indicatorWeight <- sub("SCIM_", "", indicatorWeight))
	SCIM <- within(SCIM, indicatorWeight <- gsub("__", "", indicatorWeight))
	SCIM$indicatorWeight <- as.factor(SCIM$indicatorWeight)
	SCIM$type <- NULL
	SCIM$indicator <- NULL
	SCIM$weight <- NULL
	SCIM <- within(SCIM, edgeId <- paste(fromUserId, toUserId, sep = "-"))
	SCIM$edgeId <- as.factor(SCIM$edgeId)
	SCIM$fromUserId <- NULL
	SCIM$toUserId <- NULL

	return(SCIM)
}

reshapeSCIMData <- function(SCIM) {
	# reshape so that there is one JSD column per indicator/weight pair
	SCIM <- reshape(SCIM, v.names = "distJensenShannon", timevar = "indicatorWeight",
		idvar = "edgeId", direction = "wide")
	names(SCIM) <- sub("distJensenShannon.", "", names(SCIM))
	return(SCIM)
}

processExperimentResults <- function(dir, inputFile, ignoreBaseline, experimentType) {
	sink(paste(dir, "r-anova-", experimentType, ".txt", sep = ""))

	SCIM <- reshapeSCIMData(readSCIMData(dir, inputFile, ignoreBaseline))
	indicatorWeightLevels <- tail(names(SCIM), length(names(SCIM)) - 4)

	# construct model in Wilkinson notation
	if (experimentType == "nodes") {
		factors <- c("timePeriodBase", "timePeriodLength")
	} else {
		factors <- c("timePeriodBase", "timePeriodLength", "entity")
	}
	modelStr <- paste("cbind(", paste(indicatorWeightLevels, collapse = ","),
		") ~ ", paste(factors, collapse = "*"), sep = "")
	model <- lm(as.formula(modelStr), data = SCIM)

	# Use MANOVA for within-subject effects, as it does not require sphericity. In a small
	# test case, the results appear to be consistent with univariate repeated measures
	# ANOVA.
	indicatorWeight <- factor(indicatorWeightLevels)
	manovaResult <- Anova(model, idata = data.frame(indicatorWeight),
		idesign = ~ indicatorWeight)
	manovaSummary <- summary(manovaResult, multivariate = TRUE, univariate = FALSE)
	print(manovaSummary, SSP = FALSE)

	# marginal means of factor indicatorWeight (needed for building homogenous subsets)
	margMeans <- emmeans(model, ~ rep.meas)
	print(margMeans)

	# multiple comparisons of within-subject effects (indicator/weight)
	multComp <- contrast(margMeans, method = "tukey")
	print(multComp)

	# save workspace to file
	save.image(file = paste(dir, "r-anova-", experimentType, ".RData", sep = ""))

	# plot simple effects
	for (f in factors) {
		plot <- emmip(model, as.formula(paste("~ ", f, sep = "")))
		plot <- plot + theme_linedraw() + theme(legend.position = "bottom")
		ggsave(paste(dir, experimentType, "-", f, ".pdf", sep = ""), plot = plot,
			width = 8, height = 8, units = "cm", scale = 1.5)
	}

	# plot interaction effects
	for (f1 in factors) {
		for (f2 in factors) {
			if (f1 == f2) {
				next
			}
			plot <- emmip(model, as.formula(paste(f1, " ~ ", f2, sep = "")))
			plot <- plot + theme_linedraw() + theme(legend.position = "bottom")
			ggsave(paste(dir, experimentType, "-", f2, "-", f1, ".pdf", sep = ""),
				plot = plot, width = 8, height = 8, units = "cm", scale = 1.5)
		}
	}

	sink()
}

# compare neighborhood definitions via Tukey's test (after ANOVA with single within-
# subject factor) while the other experiment parameters are kept fixed
processExperimentResultsSimple <- function(dir, inputFile, ignoreBaseline, experimentType,
		refTimeBase, refTimeLength, refEntity) {
	sink(paste(dir, "r-simple-", experimentType, ".txt", sep = ""))

	SCIM <- readSCIMData(dir, inputFile, ignoreBaseline)

	# select rows that match refTimeBase, refTimeLength, refEntity
	print(paste(refTimeBase, refTimeLength, refEntity))
	SCIM <- droplevels(SCIM[which(SCIM$timePeriodBase == refTimeBase &
		SCIM$timePeriodLength == refTimeLength & SCIM$entity == refEntity),])

	SCIM <- reshapeSCIMData(SCIM)
	print(paste("before elimination of duplicate columns:", length(names(SCIM)) - 4))
	SCIM <- SCIM[!duplicated(lapply(SCIM, summary))]
	print(paste("after:", length(names(SCIM)) - 4))
	indicatorWeightLevels <- tail(names(SCIM), length(names(SCIM)) - 4)

	# construct model in Wilkinson notation
	modelStr <- paste("cbind(", paste(indicatorWeightLevels, collapse = ","), ") ~ 1",
		sep = "")
	model <- lm(as.formula(modelStr), data = SCIM)

	# MANOVA for single within-subject factor
	indicatorWeight <- factor(indicatorWeightLevels)
	manovaResult <- Anova(model, idata = data.frame(indicatorWeight),
		idesign = ~ indicatorWeight)
	manovaSummary <- summary(manovaResult, multivariate = TRUE, univariate = FALSE)
	print(manovaSummary, SSP = FALSE)

	# marginal means of factor indicatorWeight (needed for building homogenous subsets)
	margMeans <- emmeans(model, ~ rep.meas)
	print(margMeans)

	# multiple comparisons of within-subject effects (indicator/weight)
	multComp <- contrast(margMeans, method = "tukey")
	print(multComp)

	# save workspace to file
	save.image(file = paste(dir, "r-simple-", experimentType, ".RData", sep = ""))

	sink()
}

dir <- "./"
processExperimentResults(dir, "eval-anova-edges.csv", FALSE, "edges")
processExperimentResults(dir, "eval-anova-nodes.csv", FALSE, "nodes")
#processExperimentResultsSimple(dir, "eval-anova-edges.csv", FALSE, "edges",
#	"2012.05.18", 14, "EDGE_COMMUNICATION")
#processExperimentResultsSimple(dir, "eval-anova-nodes.csv", FALSE, "nodes",
#	"2012.05.18", 14, "NODE")
