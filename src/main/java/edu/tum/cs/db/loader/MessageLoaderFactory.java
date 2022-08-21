package edu.tum.cs.db.loader;

import java.util.Set;

import edu.tum.cs.db.GenericSocialMediaDao;
import edu.tum.cs.db.SocialMediaDao;
import edu.tum.cs.db.TwitterDao;

public class MessageLoaderFactory {

	public static MessageLoader createMessageLoader(SocialMediaDao dao, Set<Long> userIds)
			throws Exception {
		if (dao instanceof TwitterDao)
			return new TwitterMessageLoader((TwitterDao) dao, userIds);
		if (dao instanceof GenericSocialMediaDao)
			return new GenericMessageLoader((GenericSocialMediaDao) dao, userIds);
		throw new IllegalArgumentException("unknown DAO '" + dao.getClass().getSimpleName() + "'");
	}

}
