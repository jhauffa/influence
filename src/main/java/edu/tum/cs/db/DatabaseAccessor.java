package edu.tum.cs.db;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.sql.Connection;
import java.sql.SQLException;
import java.util.Collection;
import java.util.logging.Logger;

import javax.sql.DataSource;

import com.carrotsearch.hppc.LongIntMap;
import com.carrotsearch.hppc.LongObjectMap;
import com.mchange.v2.c3p0.DataSources;

import edu.tum.cs.util.ExperimentConfiguration;

public class DatabaseAccessor {

	protected static final Logger logger = Logger.getLogger(DatabaseAccessor.class.getName());
	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(DatabaseAccessor.class);
	protected static final int maxFetchSize = 100000;

	private final DataSource pooledDataSource;

	protected DatabaseAccessor() throws SQLException {
		String url = cfg.getLocalProperty("URL");
		DataSource unpooledDataSource = DataSources.unpooledDataSource(url);
		pooledDataSource = DataSources.pooledDataSource(unpooledDataSource);
	}

	/**
	 * @return Connection from the pool, close the connection if you don't use it anymore!
	 */
	protected Connection getConnection() throws SQLException {
		return pooledDataSource.getConnection();
	}

	public static byte[] serializeObject(Serializable object) throws IOException {
		if ((object == null) ||
			((object instanceof Collection) && (((Collection<?>) object).size() == 0)) ||
			((object instanceof LongIntMap) && (((LongIntMap) object).size() == 0)) ||
			((object instanceof LongObjectMap) && (((LongObjectMap<?>) object).size() == 0))) {
			return null;
		}

		ByteArrayOutputStream bos = new ByteArrayOutputStream();
		ObjectOutputStream oos = new ObjectOutputStream(bos);
		try {
			oos.writeObject(object);
		} finally {
			oos.close();
		}
		return bos.toByteArray();
	}

	@SuppressWarnings("unchecked")
	public static <T> T deserializeObject(byte[] byteArray, T defaultObj) throws IOException, ClassNotFoundException {
		if ((byteArray == null) || (byteArray.length == 0))
			return defaultObj;

		ByteArrayInputStream bais = new ByteArrayInputStream(byteArray);
		ObjectInputStream ins = new ObjectInputStream(bais);
		try {
			return (T) ins.readObject();
		} finally {
			ins.close();
		}
	}

}
