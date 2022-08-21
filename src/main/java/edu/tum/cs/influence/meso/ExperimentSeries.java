package edu.tum.cs.influence.meso;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.logging.Logger;

public class ExperimentSeries implements Iterable<ModelVariantExperiment> {

	private static final Logger logger = Logger.getLogger(ExperimentSeries.class.getName());

	private final Collection<ModelVariantExperiment> experiments;

	public ExperimentSeries() {
		experiments = new ArrayList<ModelVariantExperiment>();
	}

	public ExperimentSeries(Collection<ModelVariantExperiment> experiments) {
		this.experiments = experiments;
	}

	public void addExperiment(ModelVariantExperiment experiment) {
		experiments.add(experiment);
	}

	@Override
	public Iterator<ModelVariantExperiment> iterator() {
		return experiments.iterator();
	}

	public void writePredictionResultsCsv(File f) throws IOException {
		PrintWriter w = new PrintWriter(new FileOutputStream(f));
		try {
			ModelVariantExperiment.writePredictionCsvHeader(w);
			for (ModelVariantExperiment experiment : experiments)
				experiment.writePredictionCsv(w);
		} finally {
			w.close();
		}
	}

	public void writeParameterFittingResultsCsv(File f) throws IOException {
		PrintWriter w = new PrintWriter(new FileOutputStream(f));
		try {
			ModelVariantExperiment.writeParameterFittingCsvHeader(w);
			for (ModelVariantExperiment experiment : experiments)
				experiment.writeParameterFittingCsv(w);
		} finally {
			w.close();
		}
	}

	public static ExperimentSeries readPredictionResultsCsv(File f,
			List<? extends ResultFilter<PredictionResult>> filters) throws IOException, ParseException {
		return readResultsCsv(f, true, filters, null);
	}

	public static ExperimentSeries readParameterFittingResultsCsv(File f,
			List<? extends ResultFilter<ParameterFittingResult>> filters) throws IOException, ParseException {
		return readResultsCsv(f, false, null, filters);
	}

	private static ExperimentSeries readResultsCsv(File f, boolean readPredictionResults,
			List<? extends ResultFilter<PredictionResult>> predictionFilters,
			List<? extends ResultFilter<ParameterFittingResult>> parameterFittingFilters)
					throws IOException, ParseException {
		Map<Date, ModelVariantExperiment> experiments = new HashMap<Date, ModelVariantExperiment>();
		int numLinesRead = 0;
		BufferedReader r = new BufferedReader(new FileReader(f), 16 * 1024 * 1024);
		try {
			r.readLine();	// discard column header
			String line;
			while ((line = r.readLine()) != null) {
				String[] parts = line.split(";");
				Queue<String> remainingParts = new LinkedList<String>(Arrays.asList(parts));

				ModelVariantExperiment exp;
				if (readPredictionResults)
					exp = ModelVariantExperiment.readPredictionCsv(remainingParts, predictionFilters);
				else
					exp = ModelVariantExperiment.readParameterFittingCsv(remainingParts, parameterFittingFilters);

				ModelVariantExperiment expBase = experiments.get(exp.getTimePeriodBase());
				if (expBase == null)
					experiments.put(exp.getTimePeriodBase(), exp);
				else
					expBase.merge(exp);

				if ((++numLinesRead % 1000000L) == 0)
					logger.info(numLinesRead + " lines read");
			}
		} finally {
			r.close();
		}
		return new ExperimentSeries(experiments.values());
	}

}
