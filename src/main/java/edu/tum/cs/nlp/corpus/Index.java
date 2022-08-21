package edu.tum.cs.nlp.corpus;

import java.io.Serializable;
import com.google.common.collect.BiMap;
import com.google.common.collect.HashBiMap;

public class Index<T> implements Serializable {

	private static final long serialVersionUID = -7375311438960373495L;

	protected final BiMap<T,Integer> indexMap = HashBiMap.create();
	protected int nextId = 0;
	private boolean readOnly;

	/**
	 * Adds an element to the Index and returns its index. If the element is already known, the index of this element is
	 * returned.
	 */
	public int add(T element) {
		Integer index = indexMap.get(element);
		if (index == null) {
			if (readOnly)
				throw new IllegalStateException("Set to read only, only known words are allowed. Tried to put in " +
						element);
			indexMap.put(element, nextId);
			nextId++;
			return nextId - 1;
		}
		return index;
	}

	public int getNumUniqueElements() {
		return nextId;
	}

	public int getElementId(T element) {
		Integer index = indexMap.get(element);
		if (index == null)
			return -1;
		return index;
	}

	public T getElement(int id) {
		return indexMap.inverse().get(id);
	}

	/**
	 * No elements are allowed to add, good for queries with existing bags.
	 */
	public void setReadOnly() {
		readOnly = true;
	}

}
