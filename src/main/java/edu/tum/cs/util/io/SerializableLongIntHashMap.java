package edu.tum.cs.util.io;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;

import com.carrotsearch.hppc.LongIntHashMap;

public class SerializableLongIntHashMap extends LongIntHashMap implements Serializable {

	private static final long serialVersionUID = -2677702921528051336L;

	private void readObject(ObjectInputStream stream) throws IOException, ClassNotFoundException {
		int n = stream.readInt();
		ensureCapacity(n);
		for (int i = 0; i < n; i++) {
			long key = stream.readLong();
			int value = stream.readInt();
			put(key, value);
		}
	}

	private void writeObject(ObjectOutputStream stream) throws IOException {
		stream.writeInt(assigned + (hasEmptyKey ? 1 : 0));
		for (int i = 0; i <= mask; i++) {
			if (keys[i] != 0) {
				stream.writeLong(keys[i]);
				stream.writeInt(values[i]);
			}
		}
		if (hasEmptyKey) {
			stream.writeLong(0L);
			stream.writeInt(values[mask + 1]);
		}
	}

}
