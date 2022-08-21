package edu.tum.cs.influence.meso;

import java.io.IOException;
import java.util.Date;

public interface ResultFilter<T> {

	public boolean evaluate(Date timePeriodBase, ModelVariant variant, T result);
	public void writeCsv(String path) throws IOException;

}
