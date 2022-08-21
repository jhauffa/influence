package edu.tum.cs.db;

import java.util.Collection;
import java.util.Date;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.BiMap;

import edu.tum.cs.db.entities.MessageTimestamp;
import edu.tum.cs.db.entities.SocialMediaMessage;
import edu.tum.cs.db.entities.SocialMediaUser;

public interface SocialMediaDao {

	public Set<Long> getUserIds(boolean coreOnly) throws Exception;
	public SocialMediaUser getUser(long userId) throws Exception;
	public BiMap<Long, String> getScreenNames(Set<Long> userIdSet) throws Exception;
	public List<SocialMediaMessage> getMessages(long userId, Date startDate, Date endDate) throws Exception;

	/**
	 * If type is null, any type of message is retrieved. The implementation is not required to honor minSeqLength and
	 * may return shorter sequences.
	 */
	public void buildMessageChain(MessageTimestamp.MessageType type,
			Map<MessageTimestamp, Collection<MessageTimestamp>> chain, Set<MessageTimestamp> roots,
			Set<MessageTimestamp> incompleteRoots, int minSeqLength) throws Exception;

}
