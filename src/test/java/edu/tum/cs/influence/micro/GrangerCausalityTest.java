package edu.tum.cs.influence.micro;

import static org.junit.Assert.assertFalse;

import org.junit.Test;

public class GrangerCausalityTest {

	@Test
	public void testPlaceholder() {
		assertFalse((new GrangerCausality()).isSignificantJJ(0.1));
		assertFalse((new GrangerCausality()).isSignificantJJ(1.0));
	}

}
