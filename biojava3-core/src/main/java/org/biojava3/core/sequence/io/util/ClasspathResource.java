package org.biojava3.core.sequence.io.util;

import static org.biojava3.core.sequence.io.util.IOUtils.close;
import static org.biojava3.core.sequence.io.util.IOUtils.copy;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
 
import org.biojava3.core.exceptions.ParserException;

public class ClasspathResource {

  private final String location;
  private final boolean preCache;

  public ClasspathResource(String location) {
    this(location, false);
  }

  public ClasspathResource(String location, boolean preCache) {
    this.location = location;
    this.preCache = preCache;
  }

  public InputStream getInputStream() {
    return createClasspathInputStream();
  }

  public BufferedReader getBufferedReader() {
    return new BufferedReader(new InputStreamReader(getInputStream()));
  }

  public List<String> getList() {
    return IOUtils.getList(getBufferedReader());
  }

  private InputStream createClasspathInputStream() {
    InputStream is;
    InputStream classpathIs = getClass().getClassLoader().getResourceAsStream(location);
    if(classpathIs == null) {
      throw new IllegalArgumentException("Location "+
          location+" resulted in a null InputStream");
    }
    if(preCache) {
      ByteArrayOutputStream os = new ByteArrayOutputStream();
      try {
        copy(classpathIs, os);
      } catch (IOException e) {
        throw new ParserException("Cannot copy classpath InputStream", e);
      }
      finally {
        close(classpathIs);
      }
      is = new ByteArrayInputStream(os.toByteArray());
    }
    else {
      is = classpathIs;
    }
    return is;
  }
}
