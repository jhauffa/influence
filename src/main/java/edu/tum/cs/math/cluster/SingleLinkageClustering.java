package edu.tum.cs.math.cluster;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.logging.Logger;

import edu.tum.cs.util.ExperimentConfiguration;

public abstract class SingleLinkageClustering<T> implements Serializable {

	private static final long serialVersionUID = -9011342102509214553L;
	private static final Logger logger = Logger.getLogger(SingleLinkageClustering.class.getName());
	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(SingleLinkageClustering.class);
	private static final int reportIntervalSeconds = 10;
	private static final int parallelThreshold = 1000;

	private final List<T> entities;
	private final boolean parallel;
	private int entity;
	private int[] pi;
	private double[] lambda;
	private double[] maxLambda;

	public SingleLinkageClustering(List<T> entities, boolean parallel) {
		this.entities = entities;
		this.parallel = parallel;
	}

	protected abstract double distance(T e1, T e2);

	private class ProgressMonitor extends Thread {
		@Override
		public void run() {
			while (true) {
				logger.info("entity " + entity + " of " + entities.size());
				try {
					Thread.sleep(reportIntervalSeconds * 1000);
				} catch (InterruptedException ex) {
					break;
				}
			}
		}
	}

	private class DistanceTask implements Runnable {
		private final int entity, start, end;
		private final double[] mu;

		public DistanceTask(int entity, int start, int end, double[] mu) {
			this.entity = entity;
			this.start = start;
			this.end = end;
			this.mu = mu;
		}

		@Override
		public void run() {
			T refEntity = entities.get(entity);
			for (int i = start; i < end; i++)
				mu[i] = distance(entities.get(i), refEntity);
		}
	}

	public void process() {
		int n = entities.size();
		if (n == 0)
			throw new RuntimeException("no data");

		pi = new int[n];
		lambda = new double[n];
		double[] mu = new double[n];

		int numThreads = cfg.getIntProperty(ExperimentConfiguration.PROP_NUM_THREADS,
				Runtime.getRuntime().availableProcessors());
		ExecutorService pool = null;
		Future<?>[] tasks = null;
		if (parallel) {
			pool = Executors.newFixedThreadPool(numThreads);
			tasks = new Future[numThreads];
		}

		Thread progressMonitor = new ProgressMonitor();
		progressMonitor.start();

		lambda[0] = Double.POSITIVE_INFINITY;
		for (entity = 1; entity < n; entity++) {
			pi[entity] = entity;
			lambda[entity] = Double.POSITIVE_INFINITY;

			if (parallel && (entity > parallelThreshold)) {
				int start = 0;
				int sliceSize = entity / numThreads;
				for (int i = 0; i < (numThreads - 1); i++) {
					int end = start + sliceSize;
					tasks[i] = pool.submit(new DistanceTask(entity, start, end, mu));
					start = end;
				}
				tasks[numThreads - 1] = pool.submit(new DistanceTask(entity, start, entity, mu));
				for (Future<?> task : tasks) {
					try {
						task.get();
					} catch (Exception ex) {
						throw new RuntimeException("error in worker thread", ex);
					}
				}
			} else {
				T refEntity = entities.get(entity);
				for (int j = 0; j < entity; j++)
					mu[j] = distance(entities.get(j), refEntity);
			}

			for (int j = 0; j < entity; j++) {
				if (lambda[j] >= mu[j]) {
					mu[pi[j]] = Math.min(mu[pi[j]], lambda[j]);
					lambda[j] = mu[j];
					pi[j] = entity;
				} else
					mu[pi[j]] = Math.min(mu[pi[j]], mu[j]);
			}

			for (int j = 0; j < entity; j++)
				if (lambda[j] >= lambda[pi[j]])
					pi[j] = entity;
		}

		if (pool != null)
			pool.shutdown();

		progressMonitor.interrupt();

		Set<Double> lambdaValues = new TreeSet<Double>();
		for (double l : lambda)
			lambdaValues.add(l);
		maxLambda = new double[lambdaValues.size()];
		int idx = 0;
		for (double l : lambdaValues)
			maxLambda[idx++] = l;
	}

	public int getMaxLevel() {
		return maxLambda.length - 1;
	}

	public Collection<Collection<T>> getClusters(int level) {
		Map<Integer, Collection<T>> clusters = new HashMap<Integer, Collection<T>>();
		for (int i = 0; i < lambda.length; i++) {
			// find representative entity for the cluster that contains i at the specified level
			int sigma = i;
			while (lambda[sigma] < maxLambda[level])
				sigma = pi[sigma];

			Collection<T> cluster = clusters.get(sigma);
			if (cluster == null) {
				cluster = new ArrayList<T>();
				clusters.put(sigma, cluster);
			}
			cluster.add(entities.get(i));
		}
		return clusters.values();
	}

}
