package org.biojava.bio.program.formats;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.biojava.bio.BioException;
import org.biojava.bio.EcNumber;
import org.biojava.bio.program.tagvalue.ChangeTable;
import org.biojava.bio.program.tagvalue.Formats;
import org.biojava.bio.program.tagvalue.Parser;
import org.biojava.utils.ClassTools;
import org.biojava.utils.Services;
import org.biojava.utils.lsid.LifeScienceIdentifier;
import org.biojava.utils.lsid.LifeScienceIdentifierParseException;

public class FormatTools {
  private FormatTools() {}

  private static Map LSID_2_FORMAT;

  public static final ChangeTable.Changer EC_FROM_STRING =
  new ChangeTable.Changer() {
    public Object change(Object value) {
      String sv = (String) value;
      return EcNumber.Impl.valueOf(sv);
    }
  };

  /**
   * Attempt to find aformat for a format identifer string.
   *
   * <p>The string will be resolved in the following way:
   * <ol>
   * <li>Treat the name as an LSID and search for a format class with that
   * LSID.</li>
   * <li>Load a class of that name</li>
   * <li>Load a class in the package
   * <code>org.biojava.bio.program.formats</code> with that name</li>
   * <li>Load a class in that package after replacing each '.' in the name with
   * "$" so that a search is made of inner classes.</li>
   * </ol>
   * <p>
   *
   * <p>It is not specified if the format returned is a new instance or not.</p>
   *
   * This method uses the service providor framework to find format providers.
   * If you add formats to the core biojava distribution, you must add the
   * class name to the file <code>biojava-live/resources/META-INF/services/org.biojava.bio.program.formats.Format</code>/ If you implement formats and
   * place them in your own .jar files, you should put the class name in a
   * similarly named file in your jar. This should mean that the format becomes
   * automatically registered with the system.
   *
   * @param formatName  the Stirng to use to find the format name
   * @return a Format for that name
   * @throws BioException if the format could not be resolved for some reason
   */
  public static Format getFormat(String formatName)
  throws BioException {
    // fixme: should use somethign better than BioException
    // should probaby go via jndi

    Format format = null;

    try {
      LifeScienceIdentifier lsid = LifeScienceIdentifier.valueOf(formatName);
      Map lsidResolvers = getLsid2Format();
      format = (Format) lsidResolvers.get(lsid);
    } catch (LifeScienceIdentifierParseException e) {
      // - it isn't an LSID I guess
    }

    if(format == null) {
      Class formatClass = null;

      try {
        formatClass = ClassTools.getClassLoader(Parser.class).loadClass(formatName);
      } catch (ClassNotFoundException cnfe1) {
      }

      if(formatClass == null) {
        String fn = "org.biojava.bio.program.formats." +
                    formatName;
        try {
          formatClass = ClassTools.getClassLoader(Parser.class).loadClass(fn);
        } catch (ClassNotFoundException cnfe2) {
        }
      }

      if(formatClass == null) {
        String fn = "org.biojava.bio.program.formats." +
                    formatName.replace('.', '$');
        try {
          formatClass = ClassTools.getClassLoader(Parser.class).loadClass(fn);
        } catch (ClassNotFoundException cnfe2) {
        }
      }

      if(formatClass != null) {
        try {
          format = (Format) formatClass.newInstance();
        } catch (InstantiationException e) {
          throw new BioException(
            "Could not instantiate class for name " + formatName,e );
        } catch (IllegalAccessException e) {
          throw new BioException(
            "Could not instantiate class for name " + formatName , e);
        }
      }
    }

    if(format == null) {
      throw new BioException("Could not resolve format name: " + formatName);
    }

    return format;
  }

  private static Map getLsid2Format()
  throws BioException {
    if(LSID_2_FORMAT == null) {
      try {
        LSID_2_FORMAT = new HashMap();
        ClassLoader loader = ClassTools.getClassLoader(Formats.class);

        Iterator implNames = Services.getImplementationNames(
          Format.class, loader ).iterator();

        while(implNames.hasNext()) {
          String name = (String) implNames.next();
          Class clazz = loader.loadClass(name);
          Format format = (Format) clazz.newInstance();
          LSID_2_FORMAT.put(format.getLSID(), format);
        }
      } catch (Exception e) {
        throw new BioException("Could not load service provider info for formats",e);
      }
    }

    return LSID_2_FORMAT;
  }
}
