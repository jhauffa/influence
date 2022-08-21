package edu.tum.cs.nlp.topic.model;

import java.io.Serializable;

import edu.tum.cs.nlp.corpus.Corpus;
import edu.tum.cs.nlp.corpus.ProcessedDocument;

public interface TopicModel<T extends ProcessedDocument> extends Serializable {

	public interface SamplingState extends Serializable {
		public double[][] estimatePhi();
		public int[][] getTopicAssignment();
	}

	public SamplingState train(Corpus<T> corpus);
	public double queryLogLikelihood(Corpus<T> corpus);

	public double[][] getTopicWordDistr();

}
