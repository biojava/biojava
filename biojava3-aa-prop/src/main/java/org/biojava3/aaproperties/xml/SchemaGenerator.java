package org.biojava3.aaproperties.xml;

import java.io.File;
import java.io.IOException;

import javax.xml.bind.SchemaOutputResolver;
import javax.xml.transform.Result;
import javax.xml.transform.stream.StreamResult;

public class SchemaGenerator extends SchemaOutputResolver{
	private String fileName;
	
	public SchemaGenerator(String filename){
		this.fileName = filename;
	}
	
	public Result createOutput(String namespaceUri, String suggestedFileName) throws IOException {
		File f = new File(this.fileName);
		f.createNewFile();
        return new StreamResult(f);
    }
}


