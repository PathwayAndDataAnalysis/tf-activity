package org.panda.tfactivity;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * Executes CausalPath recursively starting from the given directory, and navigating through subdirectories.
 *
 * @author Ozgun Babur
 */
public class RunRecursive
{
	public static void main(String[] args) throws IOException
	{
		run(args[0]);
	}

	private static void run(String dir) throws IOException
	{
		if (Files.isDirectory(Paths.get(dir)))
		{
			if (Files.exists(Paths.get(dir + File.separator + "parameters.txt")))
//				&&
//				!Files.exists(Paths.get(dir + File.separator + "TF-activity-results.txt")))
			{
				Main.main(new String[]{dir});
			}
			else
			{
				for (File file : new File(dir).listFiles())
				{
					if (file.isDirectory())
					{
						run(file.getPath());
					}
				}
			}
		}
	}
}
