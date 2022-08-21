package edu.tum.cs.nlp.tokenizer;

/*
 * Licensed to the Apache Software Foundation (ASF) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The ASF licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
 
/*
 TweetMotif is licensed under the Apache License 2.0: 
 http://www.apache.org/licenses/LICENSE-2.0.html
 Copyright Brendan O'Connor, Michel Krieger, and David Ahn, 2009-2010.
*/

/*
 Scala verion of TweetMotif is licensed under the Apache License 2.0: 
 http://www.apache.org/licenses/LICENSE-2.0.html
 Copyright Jason Baldridge, and David Snyder, 2011.
*/

/*
 * A direct port to Java from Scala version of 
 * Twitter tokenizer at https://bitbucket.org/jasonbaldridge/twokenize
 * Original Python version TweetMotif can be found at https://github.com/brendano/tweetmotif
 * 
 * Author: Vinh Khuc (khuc@cse.ohio-state.edu)
 * July 2011
 */ 

// Adapted to a research project at TU Munich.

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.tum.cs.util.ExperimentConfiguration;

public class TwitterTokenizer {
	
	protected static final Logger logger = Logger.getLogger(TwitterTokenizer.class.getName());
	protected static final ExperimentConfiguration cfg = new ExperimentConfiguration(TwitterTokenizer.class);
	public static final String PROP_DATA_PATH = "dataPath";
	public static final String PROP_LANGUAGE = "language";
	
	private static final Pattern Whitespace = Pattern.compile("[\\p{Z}\\s]+");
	
	private static final String punctChars = "['`\u201c\u201d\u2018\u2019\\\".?!,:;\\*/]";
	private static final String punctSeq = punctChars + "+";
	private static final String entity = "&(amp|lt|gt|quot);";
	
	//  URLs
	private static final String urlStart1  = "(https?://|www\\.)";
								// the generic expression at the end covers 2-letter country codes
	private static final String commonTLDs = "(com|co\\.uk|org|net|info|biz|gov|edu|int|mil|[a-z][a-z])";
	private static final String urlStart2  = "[A-Za-z0-9\\.-]+?\\." + commonTLDs + "(?=[/ \\W])";
	private static final String urlBody    = "[^ \\t\\r\\n<>]*?";
	private static final String urlExtraCrapBeforeEnd = "(" + punctChars + "|" + entity + ")+?";
	private static final String urlEnd     = "(\\.\\.+|[<>]|\\s|$)";
	private static final String url        = "\\b(" + urlStart1 + "|" + urlStart2 + ")" + urlBody +
											 "(?=(" + urlExtraCrapBeforeEnd + ")?" + urlEnd + ")";

	// Numeric
	private static final String timeLike   = "\\d+:\\d+";
	private static final String numNum     = "\\d+\\.\\d+";
	private static final String numberWithCommas = "(\\d+,)+?\\d{3}" + "(?=([^,]|$))";

	// includes 'Smart Quotes' (http://en.wikipedia.org/wiki/Smart_quotes)
	private static final String edgePunctChars    = "'\\\"\u201c\u201d\u2018\u2019<>\u00ab\u00bb{}" +
													"\\(\\)\\[\\]`\\^\\*\\+\\|";
	private static final String edgePunct    = "[" + edgePunctChars + "]";
	private static final String notEdgePunct = "[a-zA-Z0-9]";
								 // added punctSeq at the beginning so that ...* will also be detected
	private static final Pattern EdgePunctLeft  = Pattern.compile("(\\s|^|" + punctSeq + ")(" + edgePunct + "+)" +
																  "(" + notEdgePunct + ")");
								 // added punctSeq at the end so that *... will also be detected
	private static final Pattern EdgePunctRight = Pattern.compile("(" + notEdgePunct + ")(" + edgePunct + "+)" +
																  "(\\s|$|" + punctSeq + ")");

	// Abbreviations
	private static final String boundaryNotDot = "($|\\s|[\u201c\\\"?!,:;]|" + entity + ")";
								// e.g. " v.s. "
	private static final String aa1  = "([A-Za-z]\\.){2,}(?=" + boundaryNotDot + ")";
								// e.g. " i.e " or " R.I.P "
	private static final String aa2  = "[^A-Za-z]([A-Za-z]\\.){1,}[A-Za-z](?=" + boundaryNotDot + ")";
	private static final String standardAbbrev = "\\b([Mm]r|[Mm]rs|[Mm]s|[Dd]r|[Ss]r|[Jj]r|[Rr]ep|[Ss]en|[Ss]t)\\.";
	private static final String arbitraryAbbrev = "(" + aa1 + "|" + aa2 + "|" + standardAbbrev + ")";

	private static final String separators  = "(--+|\u2015)";
	private static final String symbols = ".*[\u266b\u00f8\u00af\u00b9\u00a8\u00b3\u00ac\u00b1\u00aa\u00d4\u00c7" +
										     "\u251c\u00a9\u2563\u00a3\u20ac\u252c\u009c\u009d\\=&]+.*";
	private static final String thingsThatSplitWords = "[^\\s\\.,]";
	private static final String embeddedApostrophe = thingsThatSplitWords + "+'" + thingsThatSplitWords + "+";

	// Emoticons
	private static final String normalEyes = "[:=]";
	private static final String wink = "[;]";
	private static final String noseArea = "(|o|O|-)"; // rather tight precision, \\S might be reasonable...
	private static final String happyMouths = "[D\\)\\]]";
	private static final String sadMouths = "[\\(\\[]";
	private static final String tongue = "[pP]";
	private static final String otherMouths = "[doO/\\\\]"; // remove forward slash if http://'s aren't cleaned

	private static final String emoticon = "(" + normalEyes + "|" + wink + ")" + noseArea +
										   "(" + tongue + "|" + otherMouths + "|" + sadMouths + "|" + happyMouths + ")";
	
	// Delimiters
	private static final Pattern Protected  = Pattern.compile(
	    "("
	    + emoticon + "|"
	    + url + "|"
	    + timeLike + "|"
	    + numNum + "|"
	    + numberWithCommas + "|"
	    + punctSeq + "|"
	    + arbitraryAbbrev + "|"
	    + separators + "|"
	    + embeddedApostrophe + ")");
	
	// Patterns for token filtering
	private static final String recipient = ".?[@]\\w+";
	private static final String nonWords = "\\W";	// for filtering out punctuation and other non-word characters
	private static final String digit = "^\\d+$";
	private static final String floatingNumber = "^\\d+[.,]\\d+$"; 

	private static final Pattern Ignored = Pattern.compile(
		"("
		+ recipient + "|"
		+ separators + "|"
		+ punctSeq + "|"
		+ symbols + "|"
		+ nonWords + "|"
		+ digit + "|"
		+ floatingNumber + ")");

	private static final Pattern Emoticon = Pattern.compile(emoticon);
	protected static final Pattern Url = Pattern.compile(url);
	private static final Pattern Topic = Pattern.compile("[#]\\w+");

	private static final String retweetIndicators = "(?:RT|via) @";

	private final Set<String> stopWords;
	private final EmojiTokenizer emojiTokenizer;

   public TwitterTokenizer() throws IOException {
	   this(new HashSet<String>());
	   String lang = cfg.getLocalProperty(PROP_LANGUAGE, "en");
	   stopWords.addAll(loadWordList(TwitterTokenizer.class.getResourceAsStream("/data/stopWords-" + lang + ".txt")));
	   File f = new File(cfg.getLocalProperty(PROP_DATA_PATH), "stopWords-" + lang + ".txt");
	   try {
		   stopWords.addAll(loadWordList(new FileInputStream(f)));
	   } catch (IOException ex) {
		   logger.log(Level.WARNING, "could not load dataset-specific stop words from '" + f.getPath() + "'", ex);
	   }
   }

   public TwitterTokenizer(Set<String> stopWords) throws IOException {
	   this.stopWords = stopWords;
	   emojiTokenizer = new EmojiTokenizer();
   }

   protected static Set<String> loadWordList(InputStream is) throws IOException {
	   Set<String> words = new HashSet<String>();
	   BufferedReader reader = new BufferedReader(new InputStreamReader(is, "UTF-8"));
	   try {
		   String line;
		   while ((line = reader.readLine()) != null)
			   if ((line.trim().length() > 0) && !line.startsWith("#"))
				   words.add(line.trim());
	   } finally {
		   reader.close();
	   }
	   return words;
   }
	
   // 'foo' => ' foo '
   private static String splitEdgePunct(String input) {
	    Matcher splitLeftMatcher  = EdgePunctLeft.matcher(input);
	    String splitLeft = splitLeftMatcher.replaceAll("$1$2 $3");
	    
	    Matcher splitRightMatcher = EdgePunctRight.matcher(splitLeft);
	    return splitRightMatcher.replaceAll("$1 $2$3");
   }
	
   // "foo   bar" => "foo bar"
   private static String squeezeWhitespace(String input) {
	   Matcher whitespaceMatcher = Whitespace.matcher(input);
	   return whitespaceMatcher.replaceAll(" ").trim();
   }
    
   protected Vector<String> filterToken(String token, boolean keepUrls) {
	   Vector<String> filteredTokens = new Vector<String>();
	   if (Url.matcher(token).find()) {
		   if (keepUrls) {
			   token = token.replaceFirst("https://t.co","http://t.co")
			   				.replaceFirst("https://bit.ly/","http://bit.ly/");
			   filteredTokens.add(token);
		   }
		   return filteredTokens;
	   } else
		   token = token.toLowerCase();
	   if (token.startsWith("http://") || token.startsWith("https://") || token.equals("http"))	// corrupted URLs
		   return filteredTokens;
	   if (stopWords.contains(token))
		   return filteredTokens;
	   if (Topic.matcher(token).matches()) {
		   String hashtag = token.substring(1);
		   if (!stopWords.contains(hashtag)) {
			   // generate two tokens for each hashtag to boost them
			   filteredTokens.add(hashtag);
			   filteredTokens.add(hashtag);
		   }
		   return filteredTokens;
	   }
	   if (Emoticon.matcher(token).matches()) {
		   filteredTokens.add(token);
		   return filteredTokens;
	   }
	   if (Ignored.matcher(token).matches())
		   return filteredTokens;
	   filteredTokens.add(token);
	   return filteredTokens;
   }
    
   // simpleTokenize should be called after using squeezeWhitespace()
   private List<String> simpleTokenize(String text, boolean keepUrls) {
	   // Do the no-brainers first
	   String excludeRT = text.replaceAll(retweetIndicators, "@");	// sorts out all "RT" tags for retweets
	   String splitPunctText = splitEdgePunct(excludeRT);
	   int textLength = splitPunctText.length();
	   
	   // Find the matches for subsequences that should be protected,
	   // e.g. URLs, 1.0, U.N.K.L.E., 12:53
	   Matcher protectedMatcher = Protected.matcher(splitPunctText);
	   
	   // The spans of the "bads" should not be split.
	   // I do this since I couldn't find an equivalent method in Java which does 
	   // the same as Protected.findAllIn(splitPunctText).matchData.toList
	   // badSpans is a vector of indices
	   Vector<Integer> badSpans = new Vector<Integer>();
	   while (protectedMatcher.find()) {
		   int start = protectedMatcher.start();
		   int end = protectedMatcher.end();
		   // if this badSpan is not empty
		   if (start != end) {
			   badSpans.add(new Integer(start));
			   badSpans.add(new Integer(end));
		   }
	   }
	   
	   // Create a list of indices to create the "goods", which can be
	   // split. We are taking "bad" spans like 
	   //     List((2,5), (8,10)) 
	   // to create 
	   //    List(0, 2, 5, 8, 10, 12)
	   // where, e.g., "12" here would be the textLength
	   Vector<Integer> indices = new Vector<Integer>();
	   // add index 0
	   indices.add(new Integer(0));
	   // add indices from badSpans
	   indices.addAll(badSpans);
	   // add index length -1 
	   indices.add(new Integer(textLength));
	   
	   Vector<Vector<String>> splitGoods = new Vector<Vector<String>>();
	   for (int i=0; i<indices.size(); i+=2) {
		   String strGood = splitPunctText.substring(indices.get(i), indices.get(i+1));
		   Vector<String> splitGood = new Vector<String>(Arrays.asList(strGood.trim().split(" ")));
		   splitGoods.add(splitGood);
	   }

	   // Storing as Vector<String> 
	   Vector<String> bads = new Vector<String>();
	   for (int i=0; i<badSpans.size(); i+=2) {
		   String strBad = splitPunctText.substring(badSpans.get(i), badSpans.get(i+1));
		   bads.add(strBad);
	   }

	   //  Re-interpolate the 'good' and 'bad' Lists, ensuring that
	   //  additional tokens from last good item get included
	   Vector<String> zippedStr = new Vector<String>();
	   if (splitGoods.size() == bads.size()) {
		   for (int i=0; i<splitGoods.size(); i++) {
			   zippedStr.addAll(splitGoods.get(i));
			   zippedStr.add(bads.get(i));
		   } 
	   } else {
		   for (int i=0; i<splitGoods.size()-1; i++) {
			   zippedStr.addAll(splitGoods.get(i));
			   zippedStr.add(bads.get(i));
		   } 		   
		   // add the last element from 'splitGoods'
		   zippedStr.addAll(splitGoods.lastElement());
	   }
	   
	   Vector<String> finalTokens = new Vector<String>();
	   for (String str: zippedStr) {
		   Vector<String> tokens = filterToken(str.trim(), keepUrls);
		   // only add non-empty tokens consisting of more than one character
		   for (String token: tokens) {
			   if (token.length() > 1)
				   finalTokens.add(token);
		   }
	   }

	   return finalTokens;
   }
   
   // the tokenize method which filters out white spaces before using simpleTokenize() 
   public List<String> tokenize(String text, boolean keepUrls) {
	   List<String> tokens = simpleTokenize(squeezeWhitespace(text), keepUrls);
	   if (!emojiTokenizer.hasData())
		   return tokens;
	   List<String> tokensAndEmoji = new ArrayList<String>(tokens.size());
	   for (String token : tokens)
		   tokensAndEmoji.addAll(emojiTokenizer.tokenize(token));
	   return tokensAndEmoji;
   }
   
}
