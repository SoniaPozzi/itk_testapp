#ifndef ITK_TESTAPPAPP_H
#define ITK_TESTAPPAPP_H

#include "MooseApp.h"

class ItkTestAppApp;

template <>
InputParameters validParams<ItkTestAppApp>();

class ItkTestAppApp : public MooseApp
{
public:
  ItkTestAppApp(InputParameters parameters);
  virtual ~ItkTestAppApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* ITK_TESTAPPAPP_H */
