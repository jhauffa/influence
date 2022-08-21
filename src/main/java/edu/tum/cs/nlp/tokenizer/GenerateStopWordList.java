package edu.tum.cs.nlp.tokenizer;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import edu.tum.cs.db.SocialMediaDao;
import edu.tum.cs.db.SocialMediaDaoFactory;
import edu.tum.cs.db.TwitterDao;
import edu.tum.cs.db.entities.SocialMediaMessage;
import edu.tum.cs.db.entities.SocialMediaUrl;
import edu.tum.cs.db.loader.MessageLoader;
import edu.tum.cs.db.loader.MessageLoaderFactory;
import edu.tum.cs.nlp.corpus.TokenizedMessage;
import edu.tum.cs.time.IntervalCalculator;
import edu.tum.cs.util.ExperimentConfiguration;

public class GenerateStopWordList {

	private static final Logger logger = Logger.getLogger(GenerateStopWordList.class.getName());
	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(TwitterTokenizer.class);

	private static final int defaultNumberOfStopWords = 200;

	private static class WordCount implements Comparable<WordCount> {
		public final String word;
		public final long count;

		public WordCount(String word, long count) {
			this.word = word;
			this.count = count;
		}

		@Override
		public int compareTo(WordCount wc) {
			return -Long.compare(count, wc.count);
		}

		@Override
		public String toString() {
			return word + ":\t" + count;
		}
	}

	private static void writeLists(Map<String, Integer> index, String outPath, String fileSuffix, int numberOfStopWords)
			throws IOException {
		int[] countN = new int[5];
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(outPath,
				"singleWords" + fileSuffix + ".txt")));
		try {
			for (String key : index.keySet()) {
				int count = index.get(key);
				if (count <= countN.length) {
					if (count == 1) {
						writer.write(key);
						writer.newLine();
					}
					countN[count - 1]++;
				}
			}
		} finally {
			writer.close();
		}
		for (int i = 0; i < countN.length; i++)
			logger.info("count" + (i + 1) + ": " + countN[i]);

		// the following code generates the actual stop word lists
		List<WordCount> wordCount = new ArrayList<WordCount>(index.size());
		for (Map.Entry<String, Integer> e : index.entrySet())
			wordCount.add(new WordCount(e.getKey(), e.getValue()));
		Collections.sort(wordCount);

		List<WordCount> words = new ArrayList<WordCount>(numberOfStopWords);
		for (int i = 0; i < numberOfStopWords; i++)
			words.add(wordCount.get(i));

		writer = new BufferedWriter(new FileWriter(new File(outPath, "stopWords" + fileSuffix + ".txt")));
		try {
			BufferedWriter writerWithWordCount = new BufferedWriter(new FileWriter(new File(outPath,
					"stopWordsWithWordCount" + fileSuffix + ".txt")));
			try {
				for (WordCount wc : words) {
					writer.write(wc.word);
					writer.newLine();

					writerWithWordCount.write(wc.toString());
					writerWithWordCount.newLine();
				}
			} finally {
				writerWithWordCount.close();
			}
		} finally {
			writer.close();
		}
	}

	private static void addToIndex(Map<String, Integer> index, String token) {
		Integer c = index.get(token);
		if (c == null)
			c = 0;
		index.put(token, c + 1);
	}

	public static void main(String[] args) throws Exception {
		int numberOfStopWords = defaultNumberOfStopWords;
		if (args.length > 0)
			numberOfStopWords = Integer.parseInt(args[0]);

		Date endDate = cfg.getDateProperty(ExperimentConfiguration.PROP_END_DATE);
		IntervalCalculator ic = new IntervalCalculator(endDate);
		Date startDate = ic.getStartDate();
		String outPath = cfg.getLocalProperty(TwitterTokenizer.PROP_DATA_PATH, ".");
		boolean tokenizeUrls = cfg.getBooleanProperty(ExperimentConfiguration.PROP_TOKENIZE_URLS, false);

		SocialMediaDao dao = SocialMediaDaoFactory.createDao();
		if (tokenizeUrls && !(dao instanceof TwitterDao)) {
			logger.warning("URL tokenization is only available for the Twitter dataset");
			return;
		}

		Set<Long> userIds = dao.getUserIds(false);
		MessageLoader loader = MessageLoaderFactory.createMessageLoader(dao, userIds);
		Map<String, Integer> wordIndex = new HashMap<String, Integer>();
		Map<String, Integer> urlIndex = new HashMap<String, Integer>();
		for (long userId : userIds) {
			List<SocialMediaMessage> origMessages = tokenizeUrls ? new ArrayList<SocialMediaMessage>() : null;
			List<TokenizedMessage> messages = loader.loadMessages(userId, startDate, endDate, origMessages);
			for (TokenizedMessage doc : messages) {
				for (String word : doc)
					addToIndex(wordIndex, word);
			}

			if (tokenizeUrls) {
				List<SocialMediaUrl> urls = ((TwitterDao) dao).getResolvedUrls(origMessages);
				for (SocialMediaUrl url : urls) {
					try {
						List<String> urlComponents = TwitterTokenizerUrlAware.tokenizeUrl(url.getResolvedUrl());
						for (String component : urlComponents)
							addToIndex(urlIndex, component);
					} catch (MalformedURLException ex) {
						logger.fine("malformed URL '" + url + "'");
					}
				}
			}
		}
		logger.info("Done with " + userIds.size() + ", word count: " + wordIndex.size() +
				", URL token count: " + urlIndex.size());

		logger.info("Starting Output");
		writeLists(wordIndex, outPath, "", numberOfStopWords);
		if (!urlIndex.isEmpty())
			writeLists(urlIndex, outPath, "Url", numberOfStopWords);

		logger.info("Finished");
	}

}
