package edu.tum.cs.time.hmm;

import java.io.File;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import edu.tum.cs.math.dist.HistogramImpl.IdentityHistogram;
import edu.tum.cs.math.dist.HistogramImpl.Log10Histogram;

public class ChainStatistics {

	public static void main(String[] args) throws Exception {
		if (args.length < 1) {
			System.err.println("usage: " + ChainStatistics.class.getSimpleName() + " chains");
			return;
		}

		int numIncompleteSequences = 0;
		int numSequences = 0;
		IdentityHistogram rawLengthHisto = new IdentityHistogram();
		IdentityHistogram zoomLengthHisto = new IdentityHistogram(10);
		DescriptiveStatistics seqLengthStats = new DescriptiveStatistics();
		Log10Histogram logDeltaHisto = new Log10Histogram();

		for (TimeSequence seq : TimeSequence.readSequences(new File(args[0]))) {
			numSequences++;
			if (seq.isIncomplete)
				numIncompleteSequences++;

			int seqLength = seq.data.size();
			rawLengthHisto.addValue(seqLength);
			zoomLengthHisto.addValue(seqLength);
			seqLengthStats.addValue(seqLength);

			for (ObservationReal obs : seq.data)
				logDeltaHisto.addValue(obs.value);
		}

		System.out.println(numSequences + " sequences, "+ (((double) numIncompleteSequences / numSequences) * 100.0) +
					"% incomplete");
		System.out.println(seqLengthStats);
		System.out.println(rawLengthHisto);
		System.out.println(zoomLengthHisto);

		System.out.println("time stamp deltas:");
		System.out.println(logDeltaHisto);
	}

}
