package edu.tum.cs.util.io;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;

import com.carrotsearch.hppc.IntObjectScatterMap;

public class SerializableIntObjectScatterMap<T> extends IntObjectScatterMap<T> implements Serializable {

	private static final long serialVersionUID = 2087933023709203656L;

	public SerializableIntObjectScatterMap() {
		super();
	}

	public SerializableIntObjectScatterMap(int initialSize, float loadFactor) {
		super(initialSize, loadFactor);
	}

	private void readObject(ObjectInputStream stream) throws IOException, ClassNotFoundException {
		int n = stream.readInt();
		ensureCapacity(n);
		for (int i = 0; i < n; i++) {
			int key = stream.readInt();
			@SuppressWarnings("unchecked")
			T value = (T) stream.readObject();
			put(key, value);
		}
	}

	private void writeObject(ObjectOutputStream stream) throws IOException {
		stream.writeInt(assigned + (hasEmptyKey ? 1 : 0));
		for (int i = 0; i <= mask; i++) {
			if (keys[i] != 0) {
				stream.writeInt(keys[i]);
				stream.writeObject(values[i]);
			}
		}
		if (hasEmptyKey) {
			stream.writeInt(0);
			stream.writeObject(values[mask + 1]);
		}
	}

}
