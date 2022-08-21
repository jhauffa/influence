package edu.tum.cs.nlp.tokenizer;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLDecoder;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;
import java.util.logging.Level;
import java.util.regex.Pattern;

import com.google.common.base.CaseFormat;
import com.google.common.base.CharMatcher;
import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import com.google.common.primitives.Ints;

import edu.tum.cs.db.entities.SocialMediaUrl;

public class TwitterTokenizerUrlAware extends TwitterTokenizer {

	private static final String dataPath = cfg.getLocalProperty(PROP_DATA_PATH);

	// the URL pattern used by the extractor
	private static final String crawlerUrl = "(?m)(http|https)\\://[a-zA-Z0-9\\-\\.]+\\.[a-zA-Z]{2,3}(:[a-zA-Z0-9]*)?" +
			"/?([a-zA-Z0-9\\-\\._\\?\\,\\'/\\\\\\+&%\\$#\\=~])*[^\\.\\,\\)\\(\\s\u2018\u2019\u201C\u201D'\"]";
	private static final Pattern CrawlerUrl = Pattern.compile(crawlerUrl);

	private final Set<String> urlStopWords;
	private final Pattern urlShortenersPattern;

	private final boolean tokenizeUrls;
	private final Map<String, String> resolvedUrls = new HashMap<String, String>();
	private final Map<String, Integer> resolvedUrlCount = new HashMap<String, Integer>();

	public TwitterTokenizerUrlAware(boolean tokenizeUrls) throws IOException, SQLException {
		super();
		this.tokenizeUrls = tokenizeUrls;

		urlStopWords = new HashSet<String>();
		try {
			urlStopWords.addAll(loadWordList(new FileInputStream(new File(dataPath, "stopWordsUrl.txt"))));
		} catch (IOException ex) {
			logger.log(Level.WARNING, "could not load web-specific stop words", ex);
		}

		Set<String> urlShorteners = loadWordList(new FileInputStream(new File(dataPath, "urlShorteners.txt")));
		urlShortenersPattern = buildUrlShortenersPattern(urlShorteners);
	}

	public TwitterTokenizerUrlAware(boolean tokenizeUrls, Set<String> stopWords, Set<String> urlStopWords,
			Set<String> urlShorteners) throws IOException, SQLException {
		super(stopWords);
		this.tokenizeUrls = tokenizeUrls;
		this.urlStopWords = urlStopWords;
		urlShortenersPattern = buildUrlShortenersPattern(urlShorteners);
	}

	private Pattern buildUrlShortenersPattern(Set<String> urlShorteners) {
		return Pattern.compile("(https?://|www\\.)(" + Joiner.on("|").join(urlShorteners) + ").*");
	}

	public void setResolvedUrls(List<SocialMediaUrl> urls) throws SQLException {
		resolvedUrls.clear();
		resolvedUrlCount.clear();
		for (SocialMediaUrl url : urls) {
			String resolvedUrl = url.getResolvedUrl();
			resolvedUrls.put(url.getOriginalUrl(), resolvedUrl);
			Integer n = resolvedUrlCount.get(resolvedUrl);
			if (n == null)
				n = 0;
			resolvedUrlCount.put(resolvedUrl, n + 1);
		}
	}

	private boolean isResolvedUrlUnique(String resolvedUrl) {
		Integer urlCount = resolvedUrlCount.get(resolvedUrl);
		if (urlCount != null)
			return (urlCount == 1);

		// should never happen
		logger.warning(resolvedUrl + " is not in resolvedUrlCount");
		return true;
	}

	@Override
	protected Vector<String> filterToken(String token, boolean keepUrls) {
		if (Url.matcher(token).find()) {
			Vector<String> filteredTokens = new Vector<String>();
			if (keepUrls) {
				String resolvedUrl = resolvedUrls.get(token);
				if (resolvedUrl != null)
					token = resolvedUrl;

				if (tokenizeUrls && !urlShortenersPattern.matcher(token).matches()) {
					try {
						List<String> urlTokens = tokenizeUrl(token);
						for (String urlToken : urlTokens) {
							if (!urlStopWords.contains(urlToken) && (urlToken.length() > 1))
								filteredTokens.add(urlToken);
						}
					} catch (MalformedURLException ex) {
						// ignore
					}
				} else if (!tokenizeUrls && CrawlerUrl.matcher(token).matches() && !isResolvedUrlUnique(token)) {
					filteredTokens.add(token);
				}
			}
			return filteredTokens;
		}
		return super.filterToken(token, keepUrls);
	}

	public static List<String> tokenizeUrl(String urlString) throws MalformedURLException {
		List<String> tokens = new ArrayList<String>();
		URL url = new URL(urlString);

		String host = url.getHost();
		if (host != null)
			tokens.add(host);
		String query = url.getQuery();
		if (query != null)
			tokens.add(query);	// TODO: should be tokenized; if there are key-value pairs, drop the key

		String path = url.getPath();
		if (path.length() > 1) {
			try {
				path = URLDecoder.decode(path, "utf-8");
			} catch (Exception ex) {	// throws IllegalArgumentException (and possibly other exceptions) for some URLs
				throw new MalformedURLException(ex.getMessage());
			}

			path = CaseFormat.UPPER_CAMEL.to(CaseFormat.LOWER_UNDERSCORE, path);
			Iterable<String> pathComponents = Splitter.on(CharMatcher.JAVA_LETTER_OR_DIGIT.negate())
					.trimResults()
					.omitEmptyStrings()
					.split(path);

			for (String pathComponent : pathComponents) {
				Integer integer = Ints.tryParse(pathComponent);
				if (integer == null)	// URLs frequently contain numbers, usually without any semantic value
					tokens.add(pathComponent);
			}
		}

		return tokens;
	}

}
