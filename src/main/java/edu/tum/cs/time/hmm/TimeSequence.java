package edu.tum.cs.time.hmm;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import be.ac.ulg.montefiore.run.jahmm.ObservationReal;

public class TimeSequence {

	private static final int incompleteSequenceMarker = -1;

	public final List<ObservationReal> data;
	public final boolean isIncomplete;

	public TimeSequence(List<ObservationReal> data, boolean isIncomplete) {
		this.data = data;
		this.isIncomplete = isIncomplete;
	}

	private static final DecimalFormat fmt = new DecimalFormat("0.###", DecimalFormatSymbols.getInstance(Locale.US));

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		if (isIncomplete)
			sb.append(incompleteSequenceMarker).append(';');
		boolean firstItem = true;
		for (ObservationReal obs : data) {
			if (firstItem)
				firstItem = false;
			else
				sb.append(';');
			sb.append(fmt.format(obs.value));
		}
		return sb.toString();
	}


	public static List<List<ObservationReal>> getObservations(List<TimeSequence> sequences) {
		List<List<ObservationReal>> observations = new ArrayList<List<ObservationReal>>(sequences.size());
		for (TimeSequence seq : sequences)
			observations.add(seq.data);
		return observations;
	}

	public static int countObservations(List<TimeSequence> sequences) {
		int numObs = 0;
		for (TimeSequence sequence : sequences)
			numObs += sequence.data.size();
		return numObs;
	}

	public static boolean[] getCompleteness(List<TimeSequence> sequences) {
		boolean[] isIncomplete = new boolean[sequences.size()];
		int idx = 0;
		for (TimeSequence seq : sequences)
			isIncomplete[idx++] = seq.isIncomplete;
		return isIncomplete;
	}

	public static PrintWriter prepareOutput(File f) throws IOException {
		return new PrintWriter(new GZIPOutputStream(new FileOutputStream(f)));
	}

	public static void writeSequences(File f, List<TimeSequence> sequences) throws IOException {
		PrintWriter w = prepareOutput(f);
		try {
			for (TimeSequence seq : sequences)
				w.println(seq);
		} finally {
			w.close();
		}
	}

	public static List<TimeSequence> readSequences(File f) throws IOException {
		BufferedReader r = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(f))));
		try {
			List<TimeSequence> sequences = new ArrayList<TimeSequence>();
			String line;
			int lineIdx = -1;
			while ((line = r.readLine()) != null) {
				lineIdx++;
				String[] parts = line.split(";");

				boolean isIncomplete = false;
				int startIdx = 0;
				if ((parts.length > 1) && (Double.parseDouble(parts[0]) == incompleteSequenceMarker)) {
					isIncomplete = true;
					startIdx = 1;
				}

				if (parts.length > startIdx) {
					List<ObservationReal> sequence = new ArrayList<ObservationReal>(parts.length - startIdx);
					for (int i = startIdx; i < parts.length; i++) {
						double delta = Double.parseDouble(parts[i]);
						if (delta < 0.0)
							throw new RuntimeException("invalid delta " + delta + " on line " + lineIdx);
						sequence.add(new ObservationReal(delta));
					}
					sequences.add(new TimeSequence(sequence, isIncomplete));
				}
			}
			return sequences;
		} finally {
			r.close();
		}
	}

}
