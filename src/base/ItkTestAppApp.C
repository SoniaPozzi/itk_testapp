#include "ItkTestAppApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

#include "ImageFunctionItk.h"


//#include "ImageFunctionItk.h"
template <>
InputParameters
validParams<ItkTestAppApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

ItkTestAppApp::ItkTestAppApp(InputParameters parameters) : MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  ItkTestAppApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  ItkTestAppApp::associateSyntax(_syntax, _action_factory);
}

ItkTestAppApp::~ItkTestAppApp() {}

// External entry point for dynamic application loading
extern "C" void
ItkTestAppApp__registerApps()
{
  ItkTestAppApp::registerApps();
}
void
ItkTestAppApp::registerApps()
{
  registerApp(ItkTestAppApp);
}

// External entry point for dynamic object registration
extern "C" void
ItkTestAppApp__registerObjects(Factory & factory)
{
  ItkTestAppApp::registerObjects(factory);
}
void
ItkTestAppApp::registerObjects(Factory & factory)
{
   registerFunction(ImageFunctionItk);
}

// External entry point for dynamic syntax association
extern "C" void
ItkTestAppApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  ItkTestAppApp::associateSyntax(syntax, action_factory);
}
void
ItkTestAppApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
