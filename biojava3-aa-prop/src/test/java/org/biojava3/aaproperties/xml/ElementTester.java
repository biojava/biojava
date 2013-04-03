package org.biojava3.aaproperties.xml;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;

import org.junit.Test;

public class ElementTester {
	@Test
	public void generateSchema() throws JAXBException, IOException{
		JAXBContext context = JAXBContext.newInstance(ElementTable.class);
		context.generateSchema(new SchemaGenerator("./src/main/resources/ElementMass.xsd"));
	}
	
	@Test
	public void readXml() throws JAXBException, IOException{
		ElementTable iTable = new ElementTable();
		// Get a JAXB Context for the object we created above
		JAXBContext jc = JAXBContext.newInstance(iTable.getClass());
		Unmarshaller u = jc.createUnmarshaller();
		iTable = (ElementTable)u.unmarshal(new FileInputStream("./src/main/resources/ElementMass.xml" ) );
		for(Element e:iTable.getElement()){
			System.out.println(e);
		}
	}
	
	@Test
	public void generateXml() throws JAXBException, IOException {
		List<Isotope> iList = new ArrayList<Isotope>();
		Element hydrogen = new Element("Hydrogen", "H", 1, null, 1.00794);
		iList.add(new Isotope("Hydrogen", 1, 1.00782503207));
		iList.add(new Isotope("Deuterium", 2, 2.0141017778));
		iList.add(new Isotope("Tritium", 3, 3.0160492777));
		hydrogen.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element helium = new Element("Helium", "He", 2, null, 4.002602);
		iList.add(new Isotope("Helium-3", 3, 3.0160293191));
		iList.add(new Isotope("Helium-4", 4, 4.00260325415));
		helium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element lithium = new Element("Lithium", "Li", 3, null, 6.941);
		iList.add(new Isotope("Lithium-6", 6, 6.015122795));
		iList.add(new Isotope("Lithium-7", 7, 7.01600455));
		lithium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element beryllium = new Element("Beryllium", "Be", 4, null, 9.012182);
		iList.add(new Isotope("Beryllium-9", 9, 9.0121822));
		beryllium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element boron = new Element("Boron", "B", 5, null, 10.811);
		iList.add(new Isotope("Boron-10", 10, 10.012937));
		iList.add(new Isotope("Boron-11", 11, 11.0093054));
		boron.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element carbon = new Element("Carbon", "C", 6, null, 12.0107);
		iList.add(new Isotope("Carbon-12", 12, 12.0));
		iList.add(new Isotope("Carbon-13", 13, 13.0033548378));
		iList.add(new Isotope("Carbon-14", 14, 14.003241989));
		carbon.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element nitrogen = new Element("Nitrogen", "N", 7, null, 14.0067);
		iList.add(new Isotope("Nitrogen-14", 14, 14.0030740048));
		iList.add(new Isotope("Nitrogen-15", 15, 15.0001088982));
		nitrogen.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element oxygen = new Element("Oxygen", "O", 8, null, 15.9994);
		iList.add(new Isotope("Oxygen-16", 16, 15.99491461956));
		iList.add(new Isotope("Oxygen-17", 17, 16.9991317));
		iList.add(new Isotope("Oxygen-18", 18, 17.999161));
		oxygen.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element fluorine = new Element("Fluorine", "F", 9, null, 18.9984032);
		iList.add(new Isotope("Fluorine-19", 19, 18.99840322));
		fluorine.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element neon = new Element("Neon", "Ne", 10, null, 20.1797);
		iList.add(new Isotope("Neon-20", 20, 19.9924401754));
		iList.add(new Isotope("Neon-21", 21, 20.99384668));
		iList.add(new Isotope("Neon-22", 22, 21.991385114));
		neon.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element sodium = new Element("Sodium", "Na", 11, null, 22.98976928);
		iList.add(new Isotope("Sodium-23", 23, 22.9897692809));
		sodium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element magnesium = new Element("Magnesium", "Mg", 12, null, 24.305);
		iList.add(new Isotope("Magnesium-24", 24, 23.9850417));
		iList.add(new Isotope("Magnesium-25", 25, 24.98583692));
		iList.add(new Isotope("Magnesium-26", 26, 25.982592929));
		magnesium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element aluminium = new Element("Aluminium", "Al", 13, null, 26.9815386);
		iList.add(new Isotope("Aluminium-27", 27, 26.98153863));
		aluminium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element silicon = new Element("Silicon", "Si", 14, null, 28.0855);
		iList.add(new Isotope("Silicon-28", 28, 27.9769265325));
		iList.add(new Isotope("Silicon-29", 29, 28.9764947));
		iList.add(new Isotope("Silicon-30", 30, 29.97377017));
		silicon.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element phosphorus = new Element("Phosphorus", "P", 15, null, 30.973762);
		iList.add(new Isotope("Phosphorus-31", 31, 30.97376163));
		phosphorus.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element sulfur = new Element("Sulfur", "S", 16, null, 32.065);
		iList.add(new Isotope("Sulfur-32", 32, 31.972071));
		iList.add(new Isotope("Sulfur-33", 33, 32.97145876));
		iList.add(new Isotope("Sulfur-34", 34, 33.9678669));
		iList.add(new Isotope("Sulfur-36", 36, 35.96708076));
		sulfur.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element chlorine = new Element("Chlorine", "Cl", 17, null, 35.453);
		iList.add(new Isotope("Chlorine-35", 35, 34.96885268));
		iList.add(new Isotope("Chlorine-37", 37, 36.96590259));
		chlorine.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element argon = new Element("Argon", "Ar", 18, null, 39.948);
		iList.add(new Isotope("Argon-36", 36, 35.967545106));
		iList.add(new Isotope("Argon-38", 38, 37.9627324));
		iList.add(new Isotope("Argon-40", 40, 39.9623831225));
		argon.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element potassium = new Element("Potassium", "K", 19, null, 39.0983);
		iList.add(new Isotope("Potassium-39", 39, 38.96370668));
		iList.add(new Isotope("Potassium-40", 40, 39.96399848));
		iList.add(new Isotope("Potassium-41", 41, 40.96182576));
		potassium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element calcium = new Element("Calcium", "Ca", 20, null, 40.078);
		iList.add(new Isotope("Calcium-40", 40, 39.96259098));
		iList.add(new Isotope("Calcium-42", 42, 41.95861801));
		iList.add(new Isotope("Calcium-43", 43, 42.9587666));
		iList.add(new Isotope("Calcium-44", 44, 43.9554818));
		iList.add(new Isotope("Calcium-46", 46, 45.9536926));
		iList.add(new Isotope("Calcium-48", 48, 47.952534));
		calcium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element scandium = new Element("Scandium", "Sc", 21, null, 44.955912);
		iList.add(new Isotope("Scandium-45", 45, 44.9559119));
		scandium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element titanium = new Element("Titanium", "Ti", 22, null, 47.867);
		iList.add(new Isotope("Titanium-46", 46, 45.9526316));
		iList.add(new Isotope("Titanium-47", 47, 46.9517631));
		iList.add(new Isotope("Titanium-48", 48, 47.9479463));
		iList.add(new Isotope("Titanium-49", 49, 48.94787));
		iList.add(new Isotope("Titanium-50", 50, 49.9447912));
		titanium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element vanadium = new Element("Vanadium", "V", 23, null, 50.9415);
		iList.add(new Isotope("Vanadium-50", 50, 49.9471585));
		iList.add(new Isotope("Vanadium-51", 51, 50.9439595));
		vanadium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element chromium = new Element("Chromium", "Cr", 24, null, 51.9961);
		iList.add(new Isotope("Chromium-50", 50, 49.9460442));
		iList.add(new Isotope("Chromium-52", 52, 51.9405075));
		iList.add(new Isotope("Chromium-53", 53, 52.9406494));
		iList.add(new Isotope("Chromium-54", 54, 53.9388804));
		chromium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element manganese = new Element("Manganese", "Mn", 25, null, 54.938045);
		iList.add(new Isotope("Manganese-55", 55, 54.9380451));
		manganese.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element iron = new Element("Iron", "Fe", 26, null, 55.845);
		iList.add(new Isotope("Iron-54", 54, 53.9396105));
		iList.add(new Isotope("Iron-56", 56, 55.9349375));
		iList.add(new Isotope("Iron-57", 57, 56.935394));
		iList.add(new Isotope("Iron-58", 58, 57.9332756));
		iron.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element cobalt = new Element("Cobalt", "Co", 27, null, 58.933195);
		iList.add(new Isotope("Cobalt-59", 59, 58.933195));
		cobalt.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element nickel = new Element("Nickel", "Ni", 28, null, 58.6934);
		iList.add(new Isotope("Nickel-58", 58, 57.9353429));
		iList.add(new Isotope("Nickel-60", 60, 59.9307864));
		iList.add(new Isotope("Nickel-61", 61, 60.931056));
		iList.add(new Isotope("Nickel-62", 62, 61.9283451));
		iList.add(new Isotope("Nickel-64", 64, 63.927966));
		nickel.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element copper = new Element("Copper", "Cu", 29, null, 63.546);
		iList.add(new Isotope("Copper-63", 63, 62.9295975));
		iList.add(new Isotope("Copper-65", 65, 64.9277895));
		copper.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element zinc = new Element("Zinc", "Zn", 30, null, 65.38);
		iList.add(new Isotope("Zinc-64", 64, 63.9291422));
		iList.add(new Isotope("Zinc-66", 66, 65.9260334));
		iList.add(new Isotope("Zinc-67", 67, 66.9271273));
		iList.add(new Isotope("Zinc-68", 68, 67.9248442));
		iList.add(new Isotope("Zinc-70", 70, 69.9253193));
		zinc.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element gallium = new Element("Gallium", "Ga", 31, null, 69.723);
		iList.add(new Isotope("Gallium-69", 69, 68.9255736));
		iList.add(new Isotope("Gallium-71", 71, 70.9247013));
		gallium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element germanium = new Element("Germanium", "Ge", 32, null, 72.64);
		iList.add(new Isotope("Germanium-70", 70, 69.9242474));
		iList.add(new Isotope("Germanium-72", 72, 71.9220758));
		iList.add(new Isotope("Germanium-73", 73, 72.9234589));
		iList.add(new Isotope("Germanium-74", 74, 73.9211778));
		iList.add(new Isotope("Germanium-76", 76, 75.9214026));
		germanium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element arsenic = new Element("Arsenic", "As", 33, null, 74.9216);
		iList.add(new Isotope("Arsenic-75", 75, 74.9215965));
		arsenic.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element selenium = new Element("Selenium", "Se", 34, null, 78.96);
		iList.add(new Isotope("Selenium-74", 74, 73.9224764));
		iList.add(new Isotope("Selenium-76", 76, 75.9192136));
		iList.add(new Isotope("Selenium-77", 77, 76.919914));
		iList.add(new Isotope("Selenium-78", 78, 77.9173091));
		iList.add(new Isotope("Selenium-80", 80, 79.9165213));
		iList.add(new Isotope("Selenium-82", 82, 81.9166994));
		selenium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element bromine = new Element("Bromine", "Br", 35, null, 79.904);
		iList.add(new Isotope("Bromine-79", 79, 78.9183371));
		iList.add(new Isotope("Bromine-81", 81, 80.9162906));
		bromine.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element krypton = new Element("Krypton", "Kr", 36, null, 83.798);
		iList.add(new Isotope("Krypton-78", 78, 77.9203648));
		iList.add(new Isotope("Krypton-80", 80, 79.916379));
		iList.add(new Isotope("Krypton-82", 82, 81.9134836));
		iList.add(new Isotope("Krypton-83", 83, 82.914136));
		iList.add(new Isotope("Krypton-84", 84, 83.911507));
		iList.add(new Isotope("Krypton-86", 86, 85.91061073));
		krypton.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element rubidium = new Element("Rubidium", "Rb", 37, null, 85.4678);
		iList.add(new Isotope("Rubidium-85", 85, 84.911789738));
		iList.add(new Isotope("Rubidium-87", 87, 86.909180527));
		rubidium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element strontium = new Element("Strontium", "Sr", 38, null, 87.62);
		iList.add(new Isotope("Strontium-84", 84, 83.913425));
		iList.add(new Isotope("Strontium-86", 86, 85.9092602));
		iList.add(new Isotope("Strontium-87", 87, 86.9088771));
		iList.add(new Isotope("Strontium-88", 88, 87.9056121));
		strontium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element yttrium = new Element("Yttrium", "Y", 39, null, 88.90585);
		iList.add(new Isotope("Yttrium-89", 89, 88.9058483));
		yttrium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element zirconium = new Element("Zirconium", "Zr", 40, null, 91.224);
		iList.add(new Isotope("Zirconium-90", 90, 89.9047044));
		iList.add(new Isotope("Zirconium-91", 91, 90.9056458));
		iList.add(new Isotope("Zirconium-92", 92, 91.9050408));
		iList.add(new Isotope("Zirconium-94", 94, 93.9063152));
		iList.add(new Isotope("Zirconium-96", 96, 95.9082734));
		zirconium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element niobium = new Element("Niobium", "Nb", 41, null, 92.90638);
		iList.add(new Isotope("Niobium-93", 93, 92.9063781));
		niobium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element molybdenum = new Element("Molybdenum", "Mo", 42, null, 95.96);
		iList.add(new Isotope("Molybdenum-92", 92, 91.906811));
		iList.add(new Isotope("Molybdenum-94", 94, 93.9050883));
		iList.add(new Isotope("Molybdenum-95", 95, 94.9058421));
		iList.add(new Isotope("Molybdenum-96", 96, 95.9046795));
		iList.add(new Isotope("Molybdenum-97", 97, 96.9060215));
		iList.add(new Isotope("Molybdenum-98", 98, 97.9054082));
		iList.add(new Isotope("Molybdenum-100", 100, 99.907477));
		molybdenum.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element technetium = new Element("Technetium", "Tc", 43, null, 98.0);
		iList.add(new Isotope("Technetium-97", 97, 96.906365));
		iList.add(new Isotope("Technetium-98", 98, 97.907216));
		iList.add(new Isotope("Technetium-99", 99, 98.9062547));
		technetium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element ruthenium = new Element("Ruthenium", "Ru", 44, null, 101.07);
		iList.add(new Isotope("Ruthenium-96", 96, 95.907598));
		iList.add(new Isotope("Ruthenium-98", 98, 97.905287));
		iList.add(new Isotope("Ruthenium-99", 99, 98.9059393));
		iList.add(new Isotope("Ruthenium-100", 100, 99.9042195));
		iList.add(new Isotope("Ruthenium-101", 101, 100.9055821));
		iList.add(new Isotope("Ruthenium-102", 102, 101.9043493));
		iList.add(new Isotope("Ruthenium-104", 104, 103.905433));
		ruthenium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element rhodium = new Element("Rhodium", "Rh", 45, null, 102.9055);
		iList.add(new Isotope("Rhodium-103", 103, 102.905504));
		rhodium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element palladium = new Element("Palladium", "Pd", 46, null, 106.42);
		iList.add(new Isotope("Palladium-102", 102, 101.905609));
		iList.add(new Isotope("Palladium-104", 104, 103.904036));
		iList.add(new Isotope("Palladium-105", 105, 104.905085));
		iList.add(new Isotope("Palladium-106", 106, 105.903486));
		iList.add(new Isotope("Palladium-108", 108, 107.903892));
		iList.add(new Isotope("Palladium-110", 110, 109.905153));
		palladium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element silver = new Element("Silver", "Ag", 47, null, 107.8682);
		iList.add(new Isotope("Silver-107", 107, 106.905097));
		iList.add(new Isotope("Silver-109", 109, 108.904752));
		silver.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element cadmium = new Element("Cadmium", "Cd", 48, null, 112.411);
		iList.add(new Isotope("Cadmium-106", 106, 105.906459));
		iList.add(new Isotope("Cadmium-108", 108, 107.904184));
		iList.add(new Isotope("Cadmium-110", 110, 109.9030021));
		iList.add(new Isotope("Cadmium-111", 111, 110.9041781));
		iList.add(new Isotope("Cadmium-112", 112, 111.9027578));
		iList.add(new Isotope("Cadmium-113", 113, 112.9044017));
		iList.add(new Isotope("Cadmium-114", 114, 113.9033585));
		iList.add(new Isotope("Cadmium-116", 116, 115.904756));
		cadmium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element indium = new Element("Indium", "In", 49, null, 114.818);
		iList.add(new Isotope("Indium-113", 113, 112.904058));
		iList.add(new Isotope("Indium-115", 115, 114.903878));
		indium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element tin = new Element("Tin", "Sn", 50, null, 118.71);
		iList.add(new Isotope("Tin-112", 112, 111.904818));
		iList.add(new Isotope("Tin-114", 114, 113.902779));
		iList.add(new Isotope("Tin-115", 115, 114.903342));
		iList.add(new Isotope("Tin-116", 116, 115.901741));
		iList.add(new Isotope("Tin-117", 117, 116.902952));
		iList.add(new Isotope("Tin-118", 118, 117.901603));
		iList.add(new Isotope("Tin-119", 119, 118.903308));
		iList.add(new Isotope("Tin-120", 120, 119.9021947));
		iList.add(new Isotope("Tin-122", 122, 121.903439));
		iList.add(new Isotope("Tin-124", 124, 123.9052739));
		tin.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element antimony = new Element("Antimony", "Sb", 51, null, 121.76);
		iList.add(new Isotope("Antimony-121", 121, 120.9038157));
		iList.add(new Isotope("Antimony-123", 123, 122.904214));
		antimony.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element tellurium = new Element("Tellurium", "Te", 52, null, 127.6);
		iList.add(new Isotope("Tellurium-120", 120, 119.90402));
		iList.add(new Isotope("Tellurium-122", 122, 121.9030439));
		iList.add(new Isotope("Tellurium-123", 123, 122.90427));
		iList.add(new Isotope("Tellurium-124", 124, 123.9028179));
		iList.add(new Isotope("Tellurium-125", 125, 124.9044307));
		iList.add(new Isotope("Tellurium-126", 126, 125.9033117));
		iList.add(new Isotope("Tellurium-128", 128, 127.9044631));
		iList.add(new Isotope("Tellurium-130", 130, 129.9062244));
		tellurium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element iodine = new Element("Iodine", "I", 53, null, 126.90447);
		iList.add(new Isotope("Iodine-127", 127, 126.904473));
		iodine.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element xenon = new Element("Xenon", "Xe", 54, null, 131.293);
		iList.add(new Isotope("Xenon-124", 124, 123.905893));
		iList.add(new Isotope("Xenon-126", 126, 125.904274));
		iList.add(new Isotope("Xenon-128", 128, 127.9035313));
		iList.add(new Isotope("Xenon-129", 129, 128.9047794));
		iList.add(new Isotope("Xenon-130", 130, 129.903508));
		iList.add(new Isotope("Xenon-131", 131, 130.9050824));
		iList.add(new Isotope("Xenon-132", 132, 131.9041535));
		iList.add(new Isotope("Xenon-134", 134, 133.9053945));
		iList.add(new Isotope("Xenon-136", 136, 135.907219));
		xenon.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element caesium = new Element("Caesium", "Cs", 55, null, 132.9054519);
		iList.add(new Isotope("Caesium-133", 133, 132.905451933));
		caesium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element barium = new Element("Barium", "Ba", 56, null, 137.327);
		iList.add(new Isotope("Barium-130", 130, 129.9063208));
		iList.add(new Isotope("Barium-132", 132, 131.9050613));
		iList.add(new Isotope("Barium-134", 134, 133.9045084));
		iList.add(new Isotope("Barium-135", 135, 134.9056886));
		iList.add(new Isotope("Barium-136", 136, 135.9045759));
		iList.add(new Isotope("Barium-137", 137, 136.9058274));
		iList.add(new Isotope("Barium-138", 138, 137.9052472));
		barium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element lanthanum = new Element("Lanthanum", "La", 57, null, 138.90547);
		iList.add(new Isotope("Lanthanum-138", 138, 137.907112));
		iList.add(new Isotope("Lanthanum-139", 139, 138.9063533));
		lanthanum.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element cerium = new Element("Cerium", "Ce", 58, null, 140.116);
		iList.add(new Isotope("Cerium-136", 136, 135.907172));
		iList.add(new Isotope("Cerium-138", 138, 137.905991));
		iList.add(new Isotope("Cerium-140", 140, 139.9054387));
		iList.add(new Isotope("Cerium-142", 142, 141.909244));
		cerium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element praseodymium = new Element("Praseodymium", "Pr", 59, null, 140.90765);
		iList.add(new Isotope("Praseodymium-141", 141, 140.9076528));
		praseodymium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element neodymium = new Element("Neodymium", "Nd", 60, null, 144.242);
		iList.add(new Isotope("Neodymium-142", 142, 141.9077233));
		iList.add(new Isotope("Neodymium-143", 143, 142.9098143));
		iList.add(new Isotope("Neodymium-144", 144, 143.9100873));
		iList.add(new Isotope("Neodymium-145", 145, 144.9125736));
		iList.add(new Isotope("Neodymium-146", 146, 145.9131169));
		iList.add(new Isotope("Neodymium-148", 148, 147.916893));
		iList.add(new Isotope("Neodymium-150", 150, 149.920891));
		neodymium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element promethium = new Element("Promethium", "Pm", 61, null, 145.0);
		iList.add(new Isotope("Promethium-145", 145, 144.912749));
		iList.add(new Isotope("Promethium-147", 147, 146.9151385));
		promethium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element samarium = new Element("Samarium", "Sm", 62, null, 150.36);
		iList.add(new Isotope("Samarium-144", 144, 143.911999));
		iList.add(new Isotope("Samarium-147", 147, 146.9148979));
		iList.add(new Isotope("Samarium-148", 148, 147.9148227));
		iList.add(new Isotope("Samarium-149", 149, 148.9171847));
		iList.add(new Isotope("Samarium-150", 150, 149.9172755));
		iList.add(new Isotope("Samarium-152", 152, 151.9197324));
		iList.add(new Isotope("Samarium-154", 154, 153.9222093));
		samarium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element europium = new Element("Europium", "Eu", 63, null, 151.964);
		iList.add(new Isotope("Europium-151", 151, 150.9198502));
		iList.add(new Isotope("Europium-153", 153, 152.9212303));
		europium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element gadolinium = new Element("Gadolinium", "Gd", 64, null, 157.25);
		iList.add(new Isotope("Gadolinium-152", 152, 151.919791));
		iList.add(new Isotope("Gadolinium-154", 154, 153.9208656));
		iList.add(new Isotope("Gadolinium-155", 155, 154.922622));
		iList.add(new Isotope("Gadolinium-156", 156, 155.9221227));
		iList.add(new Isotope("Gadolinium-157", 157, 156.9239601));
		iList.add(new Isotope("Gadolinium-158", 158, 157.9241039));
		iList.add(new Isotope("Gadolinium-160", 160, 159.9270541));
		gadolinium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element terbium = new Element("Terbium", "Tb", 65, null, 158.92535);
		iList.add(new Isotope("Terbium-159", 159, 158.9253468));
		terbium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element dysprosium = new Element("Dysprosium", "Dy", 66, null, 162.5);
		iList.add(new Isotope("Dysprosium-156", 156, 155.924283));
		iList.add(new Isotope("Dysprosium-158", 158, 157.924409));
		iList.add(new Isotope("Dysprosium-160", 160, 159.9251975));
		iList.add(new Isotope("Dysprosium-161", 161, 160.9269334));
		iList.add(new Isotope("Dysprosium-162", 162, 161.9267984));
		iList.add(new Isotope("Dysprosium-163", 163, 162.9287312));
		iList.add(new Isotope("Dysprosium-164", 164, 163.9291748));
		dysprosium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element holmium = new Element("Holmium", "Ho", 67, null, 164.93032);
		iList.add(new Isotope("Holmium-165", 165, 164.9303221));
		holmium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element erbium = new Element("Erbium", "Er", 68, null, 167.259);
		iList.add(new Isotope("Erbium-162", 162, 161.928778));
		iList.add(new Isotope("Erbium-164", 164, 163.9292));
		iList.add(new Isotope("Erbium-166", 166, 165.9302931));
		iList.add(new Isotope("Erbium-167", 167, 166.9320482));
		iList.add(new Isotope("Erbium-168", 168, 167.9323702));
		iList.add(new Isotope("Erbium-170", 170, 169.9354643));
		erbium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element thulium = new Element("Thulium", "Tm", 69, null, 168.93421);
		iList.add(new Isotope("Thulium-169", 169, 168.9342133));
		thulium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element ytterbium = new Element("Ytterbium", "Yb", 70, null, 173.054);
		iList.add(new Isotope("Ytterbium-168", 168, 167.933897));
		iList.add(new Isotope("Ytterbium-170", 170, 169.9347618));
		iList.add(new Isotope("Ytterbium-171", 171, 170.9363258));
		iList.add(new Isotope("Ytterbium-172", 172, 171.9363815));
		iList.add(new Isotope("Ytterbium-173", 173, 172.9382108));
		iList.add(new Isotope("Ytterbium-174", 174, 173.9388621));
		iList.add(new Isotope("Ytterbium-176", 176, 175.9425717));
		ytterbium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element lutetium = new Element("Lutetium", "Lu", 71, null, 174.9668);
		iList.add(new Isotope("Lutetium-175", 175, 174.9407718));
		iList.add(new Isotope("Lutetium-176", 176, 175.9426863));
		lutetium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element hafnium = new Element("Hafnium", "Hf", 72, null, 178.49);
		iList.add(new Isotope("Hafnium-174", 174, 173.940046));
		iList.add(new Isotope("Hafnium-176", 176, 175.9414086));
		iList.add(new Isotope("Hafnium-177", 177, 176.9432207));
		iList.add(new Isotope("Hafnium-178", 178, 177.9436988));
		iList.add(new Isotope("Hafnium-179", 179, 178.9458161));
		iList.add(new Isotope("Hafnium-180", 180, 179.94655));
		hafnium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element tantalum = new Element("Tantalum", "Ta", 73, null, 180.94788);
		iList.add(new Isotope("Tantalum-180", 180, 179.9474648));
		iList.add(new Isotope("Tantalum-181", 181, 180.9479958));
		tantalum.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element tungsten = new Element("Tungsten", "W", 74, null, 183.84);
		iList.add(new Isotope("Tungsten-180", 180, 179.946704));
		iList.add(new Isotope("Tungsten-182", 182, 181.9482042));
		iList.add(new Isotope("Tungsten-183", 183, 182.950223));
		iList.add(new Isotope("Tungsten-184", 184, 183.9509312));
		iList.add(new Isotope("Tungsten-186", 186, 185.9543641));
		tungsten.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element rhenium = new Element("Rhenium", "Re", 75, null, 186.207);
		iList.add(new Isotope("Rhenium-185", 185, 184.952955));
		iList.add(new Isotope("Rhenium-187", 187, 186.9557531));
		rhenium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element osmium = new Element("Osmium", "Os", 76, null, 190.23);
		iList.add(new Isotope("Osmium-184", 184, 183.9524891));
		iList.add(new Isotope("Osmium-186", 186, 185.9538382));
		iList.add(new Isotope("Osmium-187", 187, 186.9557505));
		iList.add(new Isotope("Osmium-188", 188, 187.9558382));
		iList.add(new Isotope("Osmium-189", 189, 188.9581475));
		iList.add(new Isotope("Osmium-190", 190, 189.958447));
		iList.add(new Isotope("Osmium-192", 192, 191.9614807));
		osmium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element iridium = new Element("Iridium", "Ir", 77, null, 192.217);
		iList.add(new Isotope("Iridium-191", 191, 190.960594));
		iList.add(new Isotope("Iridium-193", 193, 192.9629264));
		iridium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element platinum = new Element("Platinum", "Pt", 78, null, 195.084);
		iList.add(new Isotope("Platinum-190", 190, 189.959932));
		iList.add(new Isotope("Platinum-192", 192, 191.961038));
		iList.add(new Isotope("Platinum-194", 194, 193.9626803));
		iList.add(new Isotope("Platinum-195", 195, 194.9647911));
		iList.add(new Isotope("Platinum-196", 196, 195.9649515));
		iList.add(new Isotope("Platinum-198", 198, 197.967893));
		platinum.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element gold = new Element("Gold", "Au", 79, null, 196.966569);
		iList.add(new Isotope("Gold-197", 197, 196.9665687));
		gold.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element mercury = new Element("Mercury", "Hg", 80, null, 200.59);
		iList.add(new Isotope("Mercury-196", 196, 195.965833));
		iList.add(new Isotope("Mercury-198", 198, 197.966769));
		iList.add(new Isotope("Mercury-199", 199, 198.9682799));
		iList.add(new Isotope("Mercury-200", 200, 199.968326));
		iList.add(new Isotope("Mercury-201", 201, 200.9703023));
		iList.add(new Isotope("Mercury-202", 202, 201.970643));
		iList.add(new Isotope("Mercury-204", 204, 203.9734939));
		mercury.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element thallium = new Element("Thallium", "Tl", 81, null, 204.3833);
		iList.add(new Isotope("Thallium-203", 203, 202.9723442));
		iList.add(new Isotope("Thallium-205", 205, 204.9744275));
		thallium.setIsotopes(iList);

		iList = new ArrayList<Isotope>();
		Element lead = new Element("Lead", "Pb", 82, null, 207.2);
		iList.add(new Isotope("Lead-204", 204, 203.9730436));
		iList.add(new Isotope("Lead-206", 206, 205.9744653));
		iList.add(new Isotope("Lead-207", 207, 206.9758969));
		iList.add(new Isotope("Lead-208", 208, 207.9766521));
		lead.setIsotopes(iList);

		List<Element> eList = new ArrayList<Element>();
		eList.add(hydrogen);
		eList.add(helium);
		eList.add(lithium);
		eList.add(beryllium);
		eList.add(boron);
		eList.add(carbon);
		eList.add(nitrogen);
		eList.add(oxygen);
		eList.add(fluorine);
		eList.add(neon);
		eList.add(sodium);
		eList.add(magnesium);
		eList.add(aluminium);
		eList.add(silicon);
		eList.add(phosphorus);
		eList.add(sulfur);
		eList.add(chlorine);
		eList.add(argon);
		eList.add(potassium);
		eList.add(calcium);
		eList.add(scandium);
		eList.add(titanium);
		eList.add(vanadium);
		eList.add(chromium);
		eList.add(manganese);
		eList.add(iron);
		eList.add(cobalt);
		eList.add(nickel);
		eList.add(copper);
		eList.add(zinc);
		eList.add(gallium);
		eList.add(germanium);
		eList.add(arsenic);
		eList.add(selenium);
		eList.add(bromine);
		eList.add(krypton);
		eList.add(rubidium);
		eList.add(strontium);
		eList.add(yttrium);
		eList.add(zirconium);
		eList.add(niobium);
		eList.add(molybdenum);
		eList.add(technetium);
		eList.add(ruthenium);
		eList.add(rhodium);
		eList.add(palladium);
		eList.add(silver);
		eList.add(cadmium);
		eList.add(indium);
		eList.add(tin);
		eList.add(antimony);
		eList.add(tellurium);
		eList.add(iodine);
		eList.add(xenon);
		eList.add(caesium);
		eList.add(barium);
		eList.add(lanthanum);
		eList.add(cerium);
		eList.add(praseodymium);
		eList.add(neodymium);
		eList.add(promethium);
		eList.add(samarium);
		eList.add(europium);
		eList.add(gadolinium);
		eList.add(terbium);
		eList.add(dysprosium);
		eList.add(holmium);
		eList.add(erbium);
		eList.add(thulium);
		eList.add(ytterbium);
		eList.add(lutetium);
		eList.add(hafnium);
		eList.add(tantalum);
		eList.add(tungsten);
		eList.add(rhenium);
		eList.add(osmium);
		eList.add(iridium);
		eList.add(platinum);
		eList.add(gold);
		eList.add(mercury);
		eList.add(thallium);
		eList.add(lead);

		ElementTable iTable = new ElementTable(eList);
		// Get a JAXB Context for the object we created above
		JAXBContext context = JAXBContext.newInstance(iTable.getClass());
		// To convert ex to XML, I need a JAXB Marshaller
		Marshaller marshaller = context.createMarshaller();
		// Make the output pretty
		marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
		StringWriter sw = new StringWriter();
		// marshall the object to XML
		marshaller.marshal(iTable, sw);
		// print it out for this example
		System.out.println(sw.toString());
		BufferedWriter output = new BufferedWriter(new FileWriter("./src/main/resources/ElementMass.xml"));
		output.write(sw.toString());
		output.close();
	}

}
