package edu.tum.cs.influence.meso;

import java.util.Queue;

public class ModelVariant {

	public static enum ModelType {
		BASELINE_UNIFORM, BASELINE_RANDOM, BASELINE_PREVIOUS, CONST_COEFF, EMPTY_NEIGHBORHOOD, SCIM, SCIM_REVERSE,
		SCIM_SHUFFLE
	}

	private final int timePeriodLength;	// number of days
	private final EntityType entity;

	private ModelType type;
	private final IndicatorFunction indicator;	// ignored if type != SCIM
	private final WeightFunction weight;	// ignored if type != SCIM

	public ModelVariant(int timePeriodLength, EntityType entity, ModelType type, IndicatorFunction indicator,
			WeightFunction weight) {
		this.timePeriodLength = timePeriodLength;
		this.entity = entity;
		this.type = type;
		this.indicator = indicator;
		this.weight = weight;
	}

	public ModelVariant(ModelVariant other) {
		this.timePeriodLength = other.timePeriodLength;
		this.entity = other.entity;
		this.type = other.type;
		this.indicator = other.indicator;
		this.weight = other.weight;
	}

	public int getTimePeriodLength() {
		return timePeriodLength;
	}

	public EntityType getEntity() {
		return entity;
	}

	public ModelType getType() {
		return type;
	}

	public void setType(ModelType type) {
		this.type = type;
	}

	public IndicatorFunction getIndicator() {
		return indicator;
	}

	public WeightFunction getWeight() {
		return weight;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + timePeriodLength;
		result = prime * result + ((entity == null) ? 0 : entity.hashCode());
		result = prime * result + ((type == null) ? 0 : type.hashCode());
		result = prime * result	+ ((indicator == null) ? 0 : indicator.hashCode());
		result = prime * result + ((weight == null) ? 0 : weight.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object o) {
		if (o == this)
			return true;
		if (!(o instanceof ModelVariant))
			return false;
		ModelVariant other = (ModelVariant) o;
		return ((timePeriodLength == other.timePeriodLength) && (entity == other.entity) &&
				(type == other.type) && (indicator == other.indicator) && (weight == other.weight));
	}

	@Override
	public String toString() {
		return timePeriodLength + ";" + entity + ";" + type + ";" +
				(indicator != null ? indicator.name() : "") + ";" +
				(weight != null ? weight.name() : "");
	}

	public static String getCsvHeader() {
		return "timePeriodLength;entity;type;indicator;weight";
	}

	public static ModelVariant readCsv(Queue<String> parts) {
		int timePeriodLength = Integer.parseInt(parts.poll());
		EntityType entity = EntityType.valueOf(parts.poll());
		ModelType type = ModelType.valueOf(parts.poll());
		String indicatorStr = parts.poll();
		IndicatorFunction indicator = null;
		if (!indicatorStr.isEmpty())
			indicator = IndicatorFunction.valueOf(indicatorStr);
		String weightStr = parts.poll();
		WeightFunction weight = null;
		if (!weightStr.isEmpty())
			weight = WeightFunction.valueOf(weightStr);
		return new ModelVariant(timePeriodLength, entity, type, indicator, weight);
	}

}
