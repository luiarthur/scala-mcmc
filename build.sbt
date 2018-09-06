name := "mcmc"

version := "0.1.0"

scalaVersion := "2.12.6"

libraryDependencies ++= Seq(
  "org.scalatest" %% "scalatest" % "3.0.0" % "test"
)

// Github sbt source dependencies.
lazy val distribution = RootProject(uri("git://github.com/luiarthur/scala-distributions.git"))
lazy val root = Project("root", file(".")).dependsOn(distribution)

assemblyJarName in assembly := "mcmc.jar"
