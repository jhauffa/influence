package edu.tum.cs.nlp.corpus;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class DocumentCorpus extends Corpus<ProcessedDocument> {

	private static final long serialVersionUID = 2915483076731591765L;

	public DocumentCorpus(Index<String> bagOfWords, Collection<Document> rawDocuments) {
		super(processDocuments(rawDocuments, bagOfWords));
	}

	private DocumentCorpus(List<ProcessedDocument> documents) {
		super(documents);
	}

	private static List<ProcessedDocument> processDocuments(Collection<Document> rawDocuments,
			Index<String> bagOfWords) {
		List<ProcessedDocument> documents = new ArrayList<ProcessedDocument>(rawDocuments.size());
		int docIdx = 0;
		for (Document rawDocument : rawDocuments) {
			int[] wordIds = new int[rawDocument.size()];
			int wordIdx = 0;
			for (String word : rawDocument)
				wordIds[wordIdx++] = bagOfWords.add(word);
			documents.add(new ProcessedDocument(docIdx++, wordIds));
		}
		return documents;
	}

	@Override
	public List<Corpus<ProcessedDocument>> slice(int numSlices, boolean keepOriginal) {
		List<Corpus<ProcessedDocument>> slices = new ArrayList<Corpus<ProcessedDocument>>(numSlices);
		int documentsPerSlice = documents.size() / numSlices;
		int startIdx = 0;
		for (int i = 0; i < numSlices; i++) {
			int endIdx = Math.min(startIdx + documentsPerSlice, documents.size());
			if (i == (numSlices - 1))
				endIdx = documents.size();

			List<ProcessedDocument> sliceDocuments = documents.subList(startIdx, endIdx);
			if (!keepOriginal)
				sliceDocuments = new ArrayList<ProcessedDocument>(sliceDocuments);
			slices.add(new DocumentCorpus(sliceDocuments));

			startIdx = endIdx;
		}

		if (!keepOriginal)
			documents.clear();
		return slices;
	}

}
