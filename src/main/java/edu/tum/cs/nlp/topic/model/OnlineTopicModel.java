package edu.tum.cs.nlp.topic.model;

import java.io.Serializable;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import edu.tum.cs.nlp.corpus.Corpus;
import edu.tum.cs.nlp.corpus.ProcessedDocument;
import edu.tum.cs.util.arrays.DenseDouble2DArray;
import edu.tum.cs.util.arrays.Int2DArray;

public class OnlineTopicModel<S extends ProcessedDocument, T extends LDABasedGibbsSamplingModel<S>>
		implements Serializable {

	private static final long serialVersionUID = -3391691681137645721L;

	private final T baseModel;

	private final double[] weights;
	private final boolean normalizeWeights;

	private final LinkedList<T> modelHistory = new LinkedList<T>();
	private final LinkedList<T.SamplingState> modelStateHistory = new LinkedList<T.SamplingState>();

	/**
	 * Weights are specified in reverse order, i.e. weights[0] applies to the most recently trained model.
	 * The parameter 'normalizeWeights' determines the behavior when there are more weights than
	 * previously trained models: If 'normalizeWeights' is true, the weights for the actually present models are
	 * normalized to sum to one, and all other weights are treated as zero for this update.
	 */
	public OnlineTopicModel(T baseModel, double[] weights, boolean normalizeWeights) {
		this.baseModel = baseModel;
		this.weights = weights;
		this.normalizeWeights = normalizeWeights;
	}

	/**
	 * Rebuild model history from de-serialized models.
	 */
	public OnlineTopicModel(List<T> models, List<Corpus<S>> corpora, double[] weights, boolean normalizeWeights) {
		this(models.get(0), weights, normalizeWeights);
		Iterator<Corpus<S>> it = corpora.iterator();
		for (T model : models) {
			T.SamplingState state = (T.SamplingState) model.resumeTraining(it.next(), 0);
			modelHistory.addFirst(model);
			modelStateHistory.addFirst(state);
		}
	}

	public void update(Corpus<S> corpus) {
		if (modelHistory.isEmpty()) {
			baseModel.setKeepNumWordTopic(true);
			T.SamplingState state = (T.SamplingState) baseModel.train(corpus);
			modelHistory.addFirst(baseModel);
			modelStateHistory.addFirst(state);
		} else {
			// optionally normalize weights if there are more weights than previous models
			double[] curWeights = weights;
			int numModels = modelHistory.size();
			if (normalizeWeights && (numModels < weights.length)) {
				double sumActiveWeights = 0.0;
				for (int i = 0; i < numModels; i++)
					sumActiveWeights += weights[i];
				curWeights = new double[numModels];
				for (int i = 0; i < numModels; i++)
					curWeights[i] = weights[i] / sumActiveWeights;
			}

			// compute topic prior as weighted sum of previous topic-word distributions
			DenseDouble2DArray beta = new DenseDouble2DArray(baseModel.getBeta());
			for (int i = 0; i < curWeights.length; i++) {
				T.SamplingState curState = modelStateHistory.get(i);
				if (curState != null) {
					Int2DArray numWordTopic = curState.getNumWordTopic();
					for (int j = 0; j < beta.getNumRows(); j++)
						for (int k = 0; k < beta.getNumColumns(); k++)
							beta.adjust(j, k, curWeights[i] * numWordTopic.get(k, j));
				}
			}

			// train new model for corpus
			@SuppressWarnings("unchecked")	// Oracle's javac gives a compiler error if the cast is missing
			T model = (T) baseModel.createSimilar();
			model.setBeta(beta);
			T.SamplingState state = (T.SamplingState) model.train(corpus);
			modelHistory.addFirst(model);
			modelStateHistory.addFirst(state);
			if (modelStateHistory.size() > weights.length)
				modelStateHistory.removeLast();
		}
	}

	public T getCurrentModel() {
		return modelHistory.getFirst();
	}

	/**
	 * @return A list of all models. The latest model is the first model in the list.
	 */
	public List<T> getModelHistory() {
		return modelHistory;
	}

}
