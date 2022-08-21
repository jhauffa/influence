package edu.tum.cs.util.io;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.io.Serializable;
import java.util.Date;

import org.junit.Test;

import edu.tum.cs.db.DatabaseAccessor;

public class SerializableMapTest {

	private static final int[][] data = {
		{ 12, 3345934 }, { 344, 1923477 }, { 1, 59445 }, { 3487, 345 }, { 344545, 233 }, { 0, 3345930 }
	};

	private static <T extends Serializable> void checkSerialization(T map, boolean isEmpty)
			throws IOException, ClassNotFoundException {
		byte[] ser = DatabaseAccessor.serializeObject(map);
		T map2 = DatabaseAccessor.deserializeObject(ser, null);
		if (isEmpty)	// an empty map is represented as null by serializeObject
			map = null;
		assertEquals(map, map2);
	}

	@Test
	public void testLongIntMap() throws IOException, ClassNotFoundException {
		SerializableLongIntHashMap map = new SerializableLongIntHashMap();
		checkSerialization(map, true);
		for (int[] record : data)
			map.put(record[0], record[1]);
		assertEquals(data.length, map.size());
		checkSerialization(map, false);
	}

	@Test
	public void testLongObjectMap() throws IOException, ClassNotFoundException {
		SerializableLongObjectHashMap<Date> map = new SerializableLongObjectHashMap<Date>();
		checkSerialization(map, true);
		for (int[] record : data)
			map.put(record[0], new Date(record[1]));
		assertEquals(data.length, map.size());
		checkSerialization(map, false);
	}

	@Test
	public void testIntObjectScatterMap() throws IOException, ClassNotFoundException {
		SerializableIntObjectScatterMap<Date> map = new SerializableIntObjectScatterMap<Date>();
		checkSerialization(map, false);	// DatabaseAccessor does not know about this type
		for (int[] record : data)
			map.put(record[0], new Date(record[1]));
		assertEquals(data.length, map.size());
		checkSerialization(map, false);
	}

	@Test
	public void testIntIntScatterMap() throws IOException, ClassNotFoundException {
		SerializableIntIntScatterMap map = new SerializableIntIntScatterMap();
		checkSerialization(map, false);	// DatabaseAccessor does not know about this type
		for (int[] record : data)
			map.put(record[0], record[1]);
		assertEquals(data.length, map.size());
		checkSerialization(map, false);
	}

}
