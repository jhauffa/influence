package edu.tum.cs.util;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashSet;
import java.util.Properties;
import java.util.Set;

public class ExperimentConfiguration extends Properties {

	private static final long serialVersionUID = 7111809339044922775L;
	private static final SimpleDateFormat sdf = new SimpleDateFormat("yyyy.MM.dd");
	private static final String listSeparator = ",";

	// well-known global properties
	public static final String PROP_NUM_THREADS = "numThreads";
	public static final String PROP_END_DATE = "endDate";
	public static final String PROP_TIME_PERIODS = "timePeriods";
	public static final String PROP_TOKENIZE_URLS = "tokenizeURLs";

	public static final String PROP_NUM_TOPICS = "TopicModel.numTopics";
	public static final String PROP_TOPIC_MODEL_PATH = "TopicModel.modelPath";
	public static final String PROP_USE_ONLINE_TOPIC_MODEL = "TopicModel.online";

	private final String root;

	public ExperimentConfiguration(Class<?> cls) {
		try {
			String configFileName = System.getProperty("edu.tum.cs.util.config");
			if (configFileName != null)
				load(new FileReader(configFileName));
			else
				load(ExperimentConfiguration.class.getResourceAsStream("/experiment.properties"));
		} catch (IOException ex) {
			throw new RuntimeException("error reading configuration", ex);
		}
		root = cls.getSimpleName();
	}

	private String makeGlobal(String key) {
		return root + "." + key;
	}

	public String getProperty(String key, String defaultValue) {
		String value = super.getProperty(key);
		if (value == null) {
			value = defaultValue;
			if (value == null)
				throw new RuntimeException("required property '" + key + "' not specified");
		}
		return value;
	}

	public String getLocalProperty(String key, String defaultValue) {
		return getProperty(makeGlobal(key), defaultValue);
	}

	public String getProperty(String key) {
		return getProperty(key, null);
	}

	public String getLocalProperty(String key) {
		return getLocalProperty(key, null);
	}

	public int getIntProperty(String key, Integer defaultValue) {
		String rawValue = getProperty(key, (defaultValue != null) ? defaultValue.toString() : null);
		return Integer.parseInt(rawValue);
	}

	public int getLocalIntProperty(String key, Integer defaultValue) {
		return getIntProperty(makeGlobal(key), defaultValue);
	}

	public int getIntProperty(String key) {
		return getIntProperty(key, null);
	}

	public int getLocalIntProperty(String key) {
		return getLocalIntProperty(key, null);
	}

	public double getDoubleProperty(String key, Double defaultValue) {
		String rawValue = getProperty(key, (defaultValue != null) ? defaultValue.toString() : null);
		return Double.parseDouble(rawValue);
	}

	public double getLocalDoubleProperty(String key, Double defaultValue) {
		return getDoubleProperty(makeGlobal(key), defaultValue);
	}

	public double getDoubleProperty(String key) {
		return getDoubleProperty(key, null);
	}

	public double getLocalDoubleProperty(String key) {
		return getLocalDoubleProperty(key, null);
	}

	public boolean getBooleanProperty(String key, Boolean defaultValue) {
		String rawValue = getProperty(key, (defaultValue != null) ? defaultValue.toString() : null);
		return Boolean.parseBoolean(rawValue);
	}

	public boolean getLocalBooleanProperty(String key, Boolean defaultValue) {
		return getBooleanProperty(makeGlobal(key), defaultValue);
	}

	public boolean getBooleanProperty(String key) {
		return getBooleanProperty(key, null);
	}

	public boolean getLocalBooleanProperty(String key) {
		return getLocalBooleanProperty(key, null);
	}

	public Date getDateProperty(String key) {
		String rawValue = getProperty(key, null);
		try {
			return sdf.parse(rawValue);
		} catch (ParseException ex) {
			throw new RuntimeException("invalid date '" + rawValue + "' for key '" + key + "'");
		}
	}

	public Date getLocalDateProperty(String key) {
		return getDateProperty(makeGlobal(key));
	}

	public int[] getIntListProperty(String key) {
		String rawValue = getProperty(key, null);
		String[] parts = rawValue.split(listSeparator);
		int[] values = new int[parts.length];
		for (int i = 0; i < parts.length; i++)
			values[i] = Integer.parseInt(parts[i]);
		return values;
	}

	public int[] getLocalIntListProperty(String key) {
		return getIntListProperty(makeGlobal(key));
	}

	public static Set<Long> loadUserIds(String fileName, int maxUsers) throws IOException {
		Set<Long> idList = new HashSet<Long>();
		BufferedReader reader = new BufferedReader(new FileReader(fileName));
		try {
			int numUsers = 0;
			String line;
			while ((line = reader.readLine()) != null) {
				if (line.length() > 0)
					if (idList.add(Long.valueOf(line)))
						if (++numUsers >= maxUsers)
							break;
			}
		} finally {
			reader.close();
		}
		return idList;
	}

}
