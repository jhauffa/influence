package edu.tum.cs.util.io;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;

import com.carrotsearch.hppc.LongObjectHashMap;

public class SerializableLongObjectHashMap<T> extends LongObjectHashMap<T> implements Serializable {

	private static final long serialVersionUID = 5404950149361716927L;

	private void readObject(ObjectInputStream stream) throws IOException, ClassNotFoundException {
		int n = stream.readInt();
		ensureCapacity(n);
		for (int i = 0; i < n; i++) {
			long key = stream.readLong();
			@SuppressWarnings("unchecked")
			T value = (T) stream.readObject();
			put(key, value);
		}
	}

	private void writeObject(ObjectOutputStream stream) throws IOException {
		stream.writeInt(assigned + (hasEmptyKey ? 1 : 0));
		for (int i = 0; i <= mask; i++) {
			if (keys[i] != 0) {
				stream.writeLong(keys[i]);
				stream.writeObject(values[i]);
			}
		}
		if (hasEmptyKey) {
			stream.writeLong(0L);
			stream.writeObject(values[mask + 1]);
		}
	}

}
