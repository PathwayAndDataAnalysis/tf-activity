package org.panda.tfactivity;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class TwoGroupsSampleNames
{
	String[] control;
	String[] test;

	public TwoGroupsSampleNames(String filename) throws IOException
	{
		Set<String> controlSet = new HashSet<>();
		Set<String> testSet = new HashSet<>();

		Files.lines(Paths.get(filename)).map(l -> l.split("=")).forEach(t ->
		{
			t[0] = t[0].trim();

			if (t[0].equals("control-value-column")) controlSet.add(t[1].trim());
			else if (t[0].equals("test-value-column")) testSet.add(t[1].trim());
		});

		control = new ArrayList<>(controlSet).toArray(new String[controlSet.size()]);
		test = new ArrayList<>(testSet).toArray(new String[testSet.size()]);
	}

	public String[] getControl()
	{
		return control;
	}

	public String[] getTest()
	{
		return test;
	}

	public void filterOutMissingSamples(Set<String> available)
	{
		control = filterMissing(control, available);
		test = filterMissing(test, available);
	}

	private String[] filterMissing(String[] current, Set<String> from)
	{
		Set<String> set = new HashSet<>(Arrays.asList(current));
		set.retainAll(from);
		List<String> list = new ArrayList<>(set);
		Collections.sort(list);
		return list.toArray(new String[list.size()]);
	}
}
