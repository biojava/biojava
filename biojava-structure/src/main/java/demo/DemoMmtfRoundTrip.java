package demo;

import java.io.IOException;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.mmtf.MmtfStructureReader;
import org.biojava.nbio.structure.io.mmtf.MmtfStructureWriter;
import org.rcsb.mmtf.dataholders.MmtfBean;
import org.rcsb.mmtf.decoder.BeanToGet;
import org.rcsb.mmtf.decoder.GetToInflator;
import org.rcsb.mmtf.deserializers.MessagePackDeserializer;
import org.rcsb.mmtf.encoder.GetToBean;
import org.rcsb.mmtf.encoder.InflatorToGet;
import org.rcsb.mmtf.serializers.MessagePackSerializer;

public class DemoMmtfRoundTrip {

	public static Structure roundTrip(Structure structure) throws IOException {
		// Set up the transform from the inflator to the get api
		InflatorToGet inflatorToGet = new InflatorToGet();
		// Get the writer - this is what people implement
		MmtfStructureWriter mmtfStructureWriter = new MmtfStructureWriter(structure);
		// Do thte inflation
		mmtfStructureWriter.write(inflatorToGet);
		// Get the bean
		GetToBean getToBean = new GetToBean(inflatorToGet);
		// Serialize
		MessagePackSerializer messagePackSerializer = new MessagePackSerializer();
		byte[] byteArr =  messagePackSerializer.serialize(getToBean.getMmtfBean());	
		// Get the reader - this is the bit that people need to implement.
		MmtfStructureReader mmtfStructureReader = new MmtfStructureReader();
		// Set up the deserializer
		MessagePackDeserializer messagePackDeserializer = new MessagePackDeserializer();
		// Get the data
		MmtfBean mmtfBean = messagePackDeserializer.deserialize(byteArr);
		// Set up the data API
		BeanToGet beanToGet = new BeanToGet(mmtfBean);
		// Set up the inflator
		GetToInflator getToInflator = new GetToInflator();
		// Do the inflation
		getToInflator.read(beanToGet, mmtfStructureReader);
		// Get the structue
		return mmtfStructureReader.getStructure();
	}
}
