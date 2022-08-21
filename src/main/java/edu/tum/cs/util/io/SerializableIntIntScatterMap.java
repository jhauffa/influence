package edu.tum.cs.util.io;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;

import com.carrotsearch.hppc.IntIntScatterMap;

public class SerializableIntIntScatterMap extends IntIntScatterMap implements Serializable {

	private static final long serialVersionUID = 5278285369906055527L;

	public SerializableIntIntScatterMap() {
		super();
	}

	public SerializableIntIntScatterMap(int initialSize, float loadFactor) {
		super(initialSize, loadFactor);
	}

	private void readObject(ObjectInputStream stream) throws IOException, ClassNotFoundException {
		int n = stream.readInt();
		ensureCapacity(n);
		for (int i = 0; i < n; i++) {
			int key = stream.readInt();
			int value = stream.readInt();
			put(key, value);
		}
	}

	private void writeObject(ObjectOutputStream stream) throws IOException {
		stream.writeInt(assigned + (hasEmptyKey ? 1 : 0));
		for (int i = 0; i <= mask; i++) {
			if (keys[i] != 0) {
				stream.writeInt(keys[i]);
				stream.writeInt(values[i]);
			}
		}
		if (hasEmptyKey) {
			stream.writeInt(0);
			stream.writeInt(values[mask + 1]);
		}
	}

}
