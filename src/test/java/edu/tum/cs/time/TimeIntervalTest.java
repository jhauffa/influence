package edu.tum.cs.time;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import java.util.Date;

import org.junit.Test;

public class TimeIntervalTest {

	private static final long refTimeStamp = 12345678L;
	private static final Date refDate = new Date(refTimeStamp);

	@Test
	public void testEquality() {
		TimeInterval refInterval = new TimeInterval(refDate, 2);
		assertEquals(refInterval, refInterval);
		assertEquals(refInterval, new TimeInterval(new Date(refTimeStamp), 2));
		assertNotEquals(refInterval, new TimeInterval(new Date(refTimeStamp + 1), 2));
		assertNotEquals(refInterval, new TimeInterval(refDate, 3));

		TimeInterval t2 = new TimeInterval(refInterval.getEndDate(), -2);
		assertEquals(refInterval, t2);
	}

}
