package edu.tum.cs.nlp.corpus;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.simple.RandomSource;
import org.junit.Test;

import edu.tum.cs.nlp.topic.model.TopicModelTest.RandomCorpusBuilder;

public class CorpusTest {

	private static final int numDocuments = 1000;
	private static final int numPersons = 20;
	private static final int maxSlices = 5;
	private static final long seed = 2343919L;

	/** slicing should work properly even if the corpus is empty */
	private static <T extends ProcessedDocument> void verifySlicingEmptyCorpus(Corpus<T> corpus) {
		assertEquals(0, corpus.size());
		for (int i = 1; i <= maxSlices; i++) {
			List<Corpus<T>> slices = corpus.slice(i, true);
			assertEquals(i, slices.size());
			for (Corpus<T> slice : slices)
				assertEquals(0, slice.size());
		}
	}

	private static void verifyDocumentSlicing(Corpus<ProcessedDocument> corpus) {
		for (int i = 1; i <= maxSlices; i++) {
			int sliceSum = 0;
			for (Corpus<ProcessedDocument> slice : corpus.slice(i, true))
				sliceSum += slice.size();
			assertEquals(corpus.size(), sliceSum);
		}
	}

	private static void verifyMessageSlicing(Corpus<ProcessedMessage> corpus) {
		for (int i = 1; i <= maxSlices; i++) {
			Set<Long> prevAuthorRecipientPairs = new HashSet<Long>();
			int sliceSum = 0;
			for (Corpus<ProcessedMessage> slice : corpus.slice(i, true)) {
				sliceSum += slice.size();
				Set<Long> curAuthorRecipientPairs = new HashSet<Long>();
				for (ProcessedMessage message : slice) {
					for (int recipient : message.getRecipientIds())
						curAuthorRecipientPairs.add(((long) message.getSenderId() << 32) | recipient);
				}
				for (Long authorRecipient : curAuthorRecipientPairs)
					assertFalse(prevAuthorRecipientPairs.contains(authorRecipient));
				prevAuthorRecipientPairs.addAll(curAuthorRecipientPairs);
			}
			assertEquals(corpus.size(), sliceSum);
		}
	}

	private static void testDocumentCorpus(Corpus<? extends ProcessedDocument> corpus, RandomCorpusBuilder refCorpus) {
		assertTrue(refCorpus.bagOfWords.getNumUniqueElements() <= refCorpus.getNumWords());
		assertEquals(numDocuments, corpus.size());
		int n = 0;
		boolean[] indexSeen = new boolean[corpus.size()];
		for (ProcessedDocument doc : corpus) {
			n += doc.size();
			assertEquals(doc.size(), doc.getWordIds().length);
			assertTrue(doc.getIndex() < corpus.size());
			assertFalse(indexSeen[doc.getIndex()]);
			indexSeen[doc.getIndex()] = true;
		}
		assertEquals(n, corpus.countTokens());
	}

	private static void testMessageCorpus(Corpus<ProcessedMessage> corpus, RandomCorpusBuilder refCorpus) {
		testDocumentCorpus(corpus, refCorpus);
		assertTrue(refCorpus.bagOfPersons.getNumUniqueElements() <= numPersons);
		for (ProcessedMessage msg : corpus)
			assertTrue(msg.getRecipientIds().length >= 1);
	}

	@Test
	public void testCorpus() {
		UniformRandomProvider rand = RandomSource.create(RandomSource.XOR_SHIFT_1024_S, seed);
		RandomCorpusBuilder refCorpus = new RandomCorpusBuilder(rand, numDocuments, numPersons, 10, 5.0, 0.01, false);
		refCorpus.generateMessages(numDocuments, numPersons);

		Corpus<ProcessedDocument> docCorpus = new DocumentCorpus(new Index<String>(),
				Collections.<Document>emptyList());
		verifySlicingEmptyCorpus(docCorpus);
		docCorpus = refCorpus.buildDocumentCorpus();
		testDocumentCorpus(docCorpus, refCorpus);
		verifyDocumentSlicing(docCorpus);

		Corpus<ProcessedMessage> messageCorpus = new MessageCorpus<Long, Message<Long>>(new Index<String>(),
				new Index<Long>(), Collections.<Message<Long>>emptyList());
		verifySlicingEmptyCorpus(messageCorpus);
		messageCorpus = refCorpus.buildMessageCorpus();
		testMessageCorpus(messageCorpus, refCorpus);
		verifyMessageSlicing(messageCorpus);

		ProcessedDocumentReader<ProcessedMessage> reader = new ProcessedMessage.Reader();
		ExternalMemoryCorpus<ProcessedMessage> extMessageCorpus =
				new ExternalMemoryCorpus<ProcessedMessage>(messageCorpus, reader);
		testMessageCorpus(extMessageCorpus, refCorpus);

		// discards messages from original corpus, so test it at the very end
		SlicedCorpus<ProcessedMessage> slicedCorpus = new SlicedCorpus<ProcessedMessage>(messageCorpus, reader,
				maxSlices, false, false);
		testMessageCorpus(slicedCorpus, refCorpus);
	}

}
