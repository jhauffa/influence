package edu.tum.cs.util.io;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class Serializer {

	private static final Logger logger = Logger.getLogger(Serializer.class.getName());

	public static void saveObjectToFile(Object o, File f) {
		try {
			ObjectOutputStream out = new ObjectOutputStream(new GZIPOutputStream(new FileOutputStream(f)));
			try {
				out.writeObject(o);
			} finally {
				out.close();
			}
		} catch (Exception e) {
			logger.log(Level.SEVERE, "Could not serialize", e);
		}
	}

	@SuppressWarnings("unchecked")
	public static <T> T loadObjectFromFile(File f) {
		T o = null;
		try {
			ObjectInputStream in = new ObjectInputStream(new GZIPInputStream(
					new BufferedInputStream(new FileInputStream(f), 5 * 1024 * 1024)));
			try {
				o = (T) in.readObject();
			} finally {
				in.close();
			}
		} catch (Exception e) {
			logger.log(Level.SEVERE, "Could not de-serialize", e);
		}
		return o;
	}

}
