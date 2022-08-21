package edu.tum.cs.influence.meso.eval;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.stat.descriptive.moment.Mean;

import edu.tum.cs.db.SocialMediaDao;
import edu.tum.cs.db.SocialMediaDaoFactory;
import edu.tum.cs.db.entities.SocialMediaUser;
import edu.tum.cs.influence.meso.InfluenceExperimentsAnova;
import edu.tum.cs.math.dist.HistogramImpl.QuantizedProbabilityHistogram;
import edu.tum.cs.util.ExperimentConfiguration;

public class ComputeCompleteness {

	public static final String FILE_COMPLETENESS = "completeness.csv";

	private static void writeCompleteness(File f, Map<Long, Double> completeness, boolean includeIsolated)
			throws IOException {
		PrintWriter w = new PrintWriter(new FileOutputStream(f));
		try {
			for (Map.Entry<Long, Double> e : completeness.entrySet()) {
				double v = e.getValue();
				if (v == -1.0) {
					if (!includeIsolated)
						continue;
					v = 1.0;
				}
				w.println(e.getKey() + ";" + v);
			}
		} finally {
			w.close();
		}
	}

	private static void writeCompletenessHisto(File f, Map<Long, Double> completeness, boolean includeIsolated)
			throws IOException {
		PrintWriter w = new PrintWriter(new FileOutputStream(f));
		try {
			QuantizedProbabilityHistogram histo = new QuantizedProbabilityHistogram(10);
			Mean avgCompleteness = new Mean();
			for (Map.Entry<Long, Double> e : completeness.entrySet()) {
				double v = e.getValue();
				if (v == -1.0) {
					if (!includeIsolated)
						continue;
					v = 1.0;
				}
				histo.addValue(v);
				avgCompleteness.increment(v);
			}
			w.println(histo);
			w.println("avg. completeness;" + avgCompleteness.getResult() + ";;");
		} finally {
			w.close();
		}
	}

	public static void main(String[] args) throws Exception {
		ExperimentConfiguration cfg = new ExperimentConfiguration(InfluenceExperimentsAnova.class);
		String path = cfg.getLocalProperty(InfluenceExperimentsAnova.PROP_OUT_PATH, ".");

		boolean isSymmetric = false;
		if (args.length > 0)
			isSymmetric = Boolean.parseBoolean(args[0]);
		int minIncidentEdges = isSymmetric ? 2 : 1;

		SocialMediaDao dao = SocialMediaDaoFactory.createDao();
		Set<Long> userIdSet = dao.getUserIds(true);

		int n = 0;
		Map<Long, Double> completeness = new HashMap<Long, Double>(userIdSet.size());
		for (Long userId : userIdSet) {
			SocialMediaUser user = dao.getUser(userId);
			Set<Long> edges = user.getOutgoingEdges();
			int numIncidentEdges = edges.size();
			edges.retainAll(userIdSet);
			int numIncidentEdgesCore = edges.size();
			edges = user.getIncomingEdges();
			numIncidentEdges += edges.size();
			edges.retainAll(userIdSet);
			numIncidentEdgesCore += edges.size();

			numIncidentEdges = Math.max(0, numIncidentEdges - minIncidentEdges);
			numIncidentEdgesCore = Math.max(0, numIncidentEdgesCore - minIncidentEdges);

			if (numIncidentEdges > 0) {
				completeness.put(userId, (double) numIncidentEdgesCore / numIncidentEdges);
			} else {
				completeness.put(userId, -1.0);
				n++;
			}
		}
		System.err.println("# isolated = " + n);

		writeCompleteness(new File(path, FILE_COMPLETENESS), completeness, false);
		writeCompleteness(new File(path, "completeness-isolated.csv"), completeness, true);
		writeCompletenessHisto(new File(path, "completeness-histo.csv"), completeness, false);
		writeCompletenessHisto(new File(path, "completeness-isolated-histo.csv"), completeness, true);
	}

}
