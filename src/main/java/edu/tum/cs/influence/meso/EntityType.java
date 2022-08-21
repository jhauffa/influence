package edu.tum.cs.influence.meso;

import java.util.EnumSet;

public enum EntityType {

	EDGE_COMMUNICATION_EXPLICIT(0), EDGE_COMMUNICATION(1), NODE(2);

	public final int value;

	private EntityType(int value) {
		this.value = value;
	}

	public static EnumSet<EntityType> getApplicable(boolean isExplicit) {
		EnumSet<EntityType> entities = EnumSet.allOf(EntityType.class);
		if (!isExplicit)
			entities.remove(EDGE_COMMUNICATION_EXPLICIT);
		return entities;
	}

}
