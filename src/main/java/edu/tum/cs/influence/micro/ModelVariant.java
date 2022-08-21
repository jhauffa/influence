package edu.tum.cs.influence.micro;

import java.util.Queue;

public class ModelVariant {

	public static enum ModelType {
		SIMILARITY, CALIBRATED, GRANGER_SIM_TF, GRANGER_TF, GRANGER
	}

	public static enum InfluenceType {
		LEVEL, EDGE
	}

	public static enum NetworkType {
		FOLLOWING, MUTUAL_FOLLOWING, REPLY, MUTUAL_REPLY;
	}

	private final ModelType model;
	private final String subModel;
	private final InfluenceType influence;
	private final NetworkType network;
	private final boolean isFakeNetwork;

	public ModelVariant(ModelType model, InfluenceType influence, NetworkType network) {
		this(model, null, influence, network, false);
	}

	private ModelVariant(ModelType model, String subModel, InfluenceType influence, NetworkType network,
			boolean isFakeNetwork) {
		this.model = model;
		this.subModel = subModel;
		this.influence = influence;
		this.network = network;
		this.isFakeNetwork = isFakeNetwork;
	}

	public ModelVariant fake() {
		return new ModelVariant(model, subModel, influence, network, true);
	}

	public ModelVariant subModel(String id) {
		return new ModelVariant(model, id, influence, network, isFakeNetwork);
	}

	public InfluenceType getInfluence() {
		return influence;
	}

	public NetworkType getNetwork() {
		return network;
	}

	public boolean isFakeNetwork() {
		return isFakeNetwork;
	}

	public boolean isSubModel() {
		return (subModel != null);
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = prime + ((model == null) ? 0 : model.hashCode());
		result = prime * result + ((subModel == null) ? 0 : subModel.hashCode());
		result = prime * result + ((influence == null) ? 0 : influence.hashCode());
		result = prime * result + ((network == null) ? 0 : network.hashCode());
		return prime * result + (isFakeNetwork ? 1231 : 1237);
	}

	@Override
	public boolean equals(Object o) {
		if (o == this)
			return true;
		if (!(o instanceof ModelVariant))
			return false;
		ModelVariant other = (ModelVariant) o;
		return ((model == other.model) && ((subModel == other.subModel) || subModel.equals(other.subModel)) &&
				(influence == other.influence) && (network == other.network) && (isFakeNetwork == other.isFakeNetwork));
	}

	private String toString(String sep) {
		return model + sep + ((subModel != null) ? subModel : "") + sep + influence + sep + network +
				(isFakeNetwork ? "fake" : "");
	}

	@Override
	public String toString() {
		return toString(";");
	}

	public String getFilePrefix() {
		return toString("-");
	}

	public static String getCsvHeader() {
		return "model;subModel;influence;network";
	}

	public static ModelVariant readCsv(Queue<String> parts) {
		ModelType model = ModelType.valueOf(parts.poll());
		String subModel = parts.poll();
		if (subModel.isEmpty())
			subModel = null;
		InfluenceType influence = InfluenceType.valueOf(parts.poll());
		NetworkType network = NetworkType.valueOf(parts.poll());
		return new ModelVariant(model, subModel, influence, network, false);
	}

}
