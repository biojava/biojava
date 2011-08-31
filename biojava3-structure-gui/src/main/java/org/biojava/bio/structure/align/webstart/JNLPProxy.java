/*********************************************************************
* Uses reflection to access JNLP services.
*
* @see
*   <a target="_blank" 
*     href="http://croftsoft.com/library/tutorials/browser/">
*   Launching a Browser from Java</a>
*
* @version
*   2001-10-23
* @since
*   2001-08-31
* @author
*   <a href="http://croftsoft.com/">David Wallace Croft</a>
*********************************************************************/

package org.biojava.bio.structure.align.webstart;

import java.lang.reflect.*;
import java.net.*;


public final class  JNLPProxy
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
{

private static final Object  basicServiceObject
  = getBasicServiceObject ( );

@SuppressWarnings("rawtypes")
private static final Class   basicServiceClass
  = getBasicServiceClass ( );

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

public static void  main ( String [ ]  args )
  throws Exception
//////////////////////////////////////////////////////////////////////
{
  showDocument ( new URL ( args [ 0 ] ) );
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

@SuppressWarnings("unchecked")
public static boolean  showDocument ( URL  url )
//////////////////////////////////////////////////////////////////////
{
  if ( basicServiceObject == null )
  {
      System.out.println("basisServiceObject = null");
      return false;
  }

  try
  {
    Method  method = basicServiceClass.getMethod (
      "showDocument", new Class [ ] { URL.class } );

    Boolean  resultBoolean = ( Boolean )
      method.invoke ( basicServiceObject, new Object [ ] { url } );

    boolean success = resultBoolean.booleanValue ( );
    if ( ! success ) 
    System.out.println("invocation of method failed!");
    return success;
  }
  catch ( Exception  ex )
  {
    ex.printStackTrace ( );

    throw new RuntimeException ( ex.getMessage ( ) );
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

@SuppressWarnings({ "unchecked", "rawtypes" })
private static Object  getBasicServiceObject ( )
//////////////////////////////////////////////////////////////////////
{
  try
  {
    Class  serviceManagerClass
      = Class.forName ( "javax.jnlp.ServiceManager" );

    Method  lookupMethod = serviceManagerClass.getMethod ( "lookup",
      new Class [ ] { String.class } );

    return lookupMethod.invoke (
      null, new Object [ ] { "javax.jnlp.BasicService" } );
  }
  catch ( Exception  ex )
  {
    return null;
  }
}

@SuppressWarnings("rawtypes")
private static Class  getBasicServiceClass ( )
//////////////////////////////////////////////////////////////////////
{
  try
  {
    return Class.forName ( "javax.jnlp.BasicService" );
  }
  catch ( Exception  ex )
  {
    return null;
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

private  JNLPProxy ( ) { }

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
}

