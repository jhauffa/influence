package edu.tum.cs.math.dist;

import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

public abstract class Histogram<T> implements Iterable<HistogramBin<T>> {

	private final Map<Integer, Long> bins = new TreeMap<Integer, Long>();
	private long n = 0;

	protected abstract int assignBin(T value);
	protected abstract T getLowerLimit(int binIdx);
	protected abstract T getUpperLimit(int binIdx);

	public int addValue(T value) {
		int binIdx = assignBin(value);
		Long c = bins.get(binIdx);
		if (c == null)
			c = 0L;
		bins.put(binIdx, c + 1);
		n++;
		return binIdx;
	}

	public long getSampleSize() {
		return n;
	}

	public Iterator<HistogramBin<T>> iterator() {
		return new Iterator<HistogramBin<T>>() {
			private final Iterator<Map.Entry<Integer, Long>> it = bins.entrySet().iterator();

			@Override
			public boolean hasNext() {
				return it.hasNext();
			}

			@Override
			public HistogramBin<T> next() {
				Map.Entry<Integer, Long> e = it.next();
				return new HistogramBin<T>(e.getKey(), getLowerLimit(e.getKey()), getUpperLimit(e.getKey()),
						e.getValue(), (double) e.getValue() / n);
			}

			@Override
			public void remove() {
				throw new UnsupportedOperationException();
			}
		};
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("L < x;x <= U;count;prob.\n");
		for (HistogramBin<T> bin : this) {
			sb.append(bin.lowerLimit).append(';').append(bin.upperLimit).append(';').append(bin.count).append(';');
			sb.append(bin.p).append('\n');
		}
		sb.append(n).append(" values;;;\n");
		return sb.toString();
	}

}
