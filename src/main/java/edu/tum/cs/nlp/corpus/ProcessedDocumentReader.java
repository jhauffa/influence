package edu.tum.cs.nlp.corpus;

import java.io.IOException;
import java.io.ObjectInputStream;

public interface ProcessedDocumentReader<P extends ProcessedDocument> {

	public P read(ObjectInputStream is) throws IOException;

}
