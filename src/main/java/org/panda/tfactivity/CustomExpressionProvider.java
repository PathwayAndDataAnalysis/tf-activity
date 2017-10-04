package org.panda.tfactivity;

import org.panda.resource.tcga.CustomExpressionReader;

import java.io.FileNotFoundException;
import java.util.HashMap;

/**
 * @author Ozgun Babur
 */
public class CustomExpressionProvider extends TCGAExpressionProvider
{
	public CustomExpressionProvider(String file) throws FileNotFoundException
	{
		er = new CustomExpressionReader(file);
		cache = new HashMap<>();
	}
}
