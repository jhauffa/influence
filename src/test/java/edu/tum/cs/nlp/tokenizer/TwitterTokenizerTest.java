package edu.tum.cs.nlp.tokenizer;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

import org.junit.Test;

public class TwitterTokenizerTest {

	private static final String[] testText = {
		"abc ABC ABc abC",
		"abc abc  abc   abc\tabc\nabc\r\nabc abc\u00a0abc",
		"123 abc abc 1 12 12.34 5,678",
		"don't won't you're foo'bar",
		"he said \"something or someone\" foo",
		"He said, she said.",
		"I am a cat",
		"\u201ccogito\u201d ergo \u2018sum\u2019",
		"http://www.abc.com http://abc.com www.abc.com www.def.co.jp www.ghi.gov http:// http **www.example.com**",
		"...what the!?",
		"...*so*...",
		"Mrs. v.s. i.e S.M.P 13:37",
		"abc-def abc--def abc---def",
		"\u266b \u266b\u266b \u266bxyz xyz\u266bxyz drum&bass",
		":-) :) :P :p =D",
		"&gt; &lt;3 drum&amp;bass &quot;lol&quot; &foo; abc&foo;def",
		"RT @BillGates lol @SteveJobs blabla via @SteveBallmer",
		"#YourSong",
		"some 60\u009d text \u009c20",
		"\ud83e\udd3c\ud83e\udd3c\u26f9\ufe0f\u200d\u2640\ufe0f\ud83e\udd3c"
	};

	private static final String[][] testRefTokens = {
		{ "abc", "abc", "abc", "abc" },	// all letters are lower-cased
		{ "abc", "abc", "abc", "abc", "abc", "abc", "abc", "abc", "abc" },	// all kinds of whitespace are stripped
		{ "abc", "abc" },	// numbers are removed
		{ "don't", "won't", "you're", "foo'bar" },	// no special treatment of common contractions
		{ "he", "said", "something", "or", "someone", "foo" },	// quotation marks are removed
		{ "he", "said", "she", "said" },	// punctuation is removed
		{ "am", "cat" },	// tokens of length 1 are removed
		{ "cogito", "ergo", "sum" },	// fancy quotes are removed
		{ },	// URLs are removed
		{ "what", "the" },	// punctuation sequences at word boundaries are removed
		{ "so" },
		{ "mrs.", "v.s.", "i.e", "s.m.p", "13:37" },	// abbreviations and time values are not split up
		{ "abc-def", "abc", "def", "abc", "def" },	// punctuation sequences are removed
		{ },	// any token that contains a symbol is removed
		{ ":-)", ":)", ":p", ":p", "=d" },	// emoticons are preserved, but letters are lower-cased
		{ "bass", ";d", "ef" },	// special treatment of HTML entities disabled -> inconsistent results
		{ "lol", "blabla" },	// @-mentions are removed, "RT/via @" is normalized to simple @-mention
		{ "yoursong", "yoursong" },	// hashtags are duplicated

		/* Special case: Some misguided Twitter bot puts \u009c and \u009d in front of and after numbers (presumably
		   instead of a currency sign). Apparently, in earlier configurations of our software stack, this character got
		   discarded at some point, while now it is handled correctly. For compatibility with earlier bags-of-words, we
		   have to ensure it gets removed early, so that the remaining numeral is removed as well. */
		{ "some", "text" },

		/* Emoji are preserved; even if not separated by whitespace, each emoji is a separate token. Complex emoji (skin
		   tone and gender variants, etc.) are handled correctly. */
		{ "\ud83e\udd3c", "\ud83e\udd3c", "\u26f9\ufe0f\u200d\u2640\ufe0f", "\ud83e\udd3c" }
	};

	@Test
	public void testTokenization() throws Exception {
		TwitterTokenizer tokenizer = new TwitterTokenizer(Collections.<String>emptySet());
		int idx = 0;
		for (String text : testText) {
			List<String> tokens = tokenizer.tokenize(text, false);
			assertEquals(Arrays.asList(testRefTokens[idx++]), tokens);
		}
	}

	@Test
	public void testStopWords() throws Exception {
		TwitterTokenizer tokenizer = new TwitterTokenizer(new HashSet<String>(Arrays.asList("bla")));
		List<String> tokens = tokenizer.tokenize("Test bla bla blub", false);
		assertEquals(Arrays.asList(new String[] { "test", "blub" }), tokens);
		tokens = tokenizer.tokenize("#bla", false);	// stopwords also affect hashtags
		assertEquals(Collections.emptyList(), tokens);
	}

	@Test
	public void testKeepUrls() throws Exception {
		TwitterTokenizer tokenizer = new TwitterTokenizer(Collections.<String>emptySet());
		List<String> tokens = tokenizer.tokenize("http://www.example.com", false);
		assertEquals(Collections.emptyList(), tokens);
		tokens = tokenizer.tokenize("http://www.example.com", true);
		assertEquals(Arrays.asList(new String[] { "http://www.example.com" }), tokens);

		// URL prefixes of frequently used shorteners are normalized; within URLs, case is preserved
		tokens = tokenizer.tokenize("https://t.co/AaAa", true);
		assertEquals(Arrays.asList(new String[] { "http://t.co/AaAa" }), tokens);
		tokens = tokenizer.tokenize("https://bit.ly/BbBb", true);
		assertEquals(Arrays.asList(new String[] { "http://bit.ly/BbBb" }), tokens);
	}

}
