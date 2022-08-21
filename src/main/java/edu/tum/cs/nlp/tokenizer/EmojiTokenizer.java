package edu.tum.cs.nlp.tokenizer;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class EmojiTokenizer {

	private static class CharTrie {
		private static final char endToken = Character.MAX_VALUE;

		private Map<Character, CharTrie> children = null;

		public CharTrie add(char v) {
			if (children == null)
				children = new HashMap<Character, CharTrie>();
			CharTrie node = children.get(v);
			if (node == null) {
				node = new CharTrie();
				children.put(v, node);
			}
			return node;
		}

		public boolean isEmpty() {
			return (children == null);
		}

		public CharTrie add(int v) {
			CharTrie node = this;
			for (char c : Character.toChars(v))
				node = node.add(c);
			return node;
		}

		public void addEnd() {
			add(endToken);
		}

		public int match(CharSequence seq, int offset) {
			CharTrie node = this;
			int idx = offset;
			int maxEndIdx = idx;
			while ((idx < seq.length()) && (node.children != null)) {
				node = node.children.get(seq.charAt(idx));
				if (node == null)
					break;
				idx++;
				if (node.children.get(endToken) != null)
					maxEndIdx = idx;
			}
			return maxEndIdx;
		}
	}

	private CharTrie loadUnicodeSequences(InputStream is) throws IOException {
		CharTrie trie = new CharTrie();
		BufferedReader r = new BufferedReader(new InputStreamReader(is, "UTF-8"));
		try {
			String line;
			while ((line = r.readLine()) != null) {
				if (line.isEmpty() || line.startsWith("#"))
					continue;
				int commentStart = line.indexOf(';');
				if (commentStart >= 0)
					line = line.substring(0, commentStart);
				String[] parts = line.split(" ");
				if (parts.length == 1) {
					int ellipsisIdx = parts[0].indexOf("..");
					if (ellipsisIdx > 0) {
						int firstCodePoint = Integer.parseInt(parts[0].substring(0, ellipsisIdx), 16);
						int lastCodePoint = Integer.parseInt(parts[0].substring(ellipsisIdx + 2), 16);
						for (int i = firstCodePoint; i <= lastCodePoint; i++)
							trie.add(i).addEnd();
					} else {
						trie.add(Integer.parseInt(parts[0], 16)).addEnd();
					}
				} else {
					CharTrie curNode = trie;
					for (String part : parts)
						curNode = curNode.add(Integer.parseInt(part, 16));
					curNode.addEnd();
				}
			}
		} finally {
			r.close();
		}
		return trie;
	}

	private final CharTrie emojiTrie;

	public EmojiTokenizer() throws IOException {
		emojiTrie = loadUnicodeSequences(EmojiTokenizer.class.getResourceAsStream("/data/emoji.txt"));
	}

	public boolean hasData() {
		return !emojiTrie.isEmpty();
	}

	public List<String> tokenize(String text) {
		List<String> tokens = new ArrayList<String>();
		int tokenStart = 0, pos = 0, tokenEnd;
		while (pos < text.length()) {
			tokenEnd = emojiTrie.match(text, pos);
			if (tokenEnd == pos) {
				pos++;
			} else {
				if (tokenStart < pos)
					tokens.add(text.substring(tokenStart, pos));
				tokens.add(text.substring(pos, tokenEnd));
				tokenStart = pos = tokenEnd;
			}
		}
		if (tokenStart < pos)
			tokens.add(text.substring(tokenStart, pos));
		return tokens;
	}

}
