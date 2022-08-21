package edu.tum.cs.influence.micro;

public class WeightedEdge<T> {

	public final int id;
	public final T weight;

	public WeightedEdge(int id, T weight) {
		this.id = id;
		this.weight = weight;
	}

	@Override
	public int hashCode() {
		return 31 + id;
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof WeightedEdge))
			return false;
		return (id == ((WeightedEdge<?>) obj).id);
	}

}
