package edu.tum.cs.db;

import edu.tum.cs.util.ExperimentConfiguration;

public class SocialMediaDaoFactory {

	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(SocialMediaDaoFactory.class);

	public static SocialMediaDao createDao() throws Exception {
		String dataSource = cfg.getLocalProperty("dataSource");
		if (dataSource.equals("twitter"))
			return new TwitterDao();
		return new GenericSocialMediaDao(dataSource);
	}

}
