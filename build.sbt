name := "mcmc"

version := "0.1.0"

scalaVersion := "2.12.6"

libraryDependencies ++= Seq(
  "org.apache.commons" % "commons-math3" % "3.6.1",
  //"org.ddahl" %% "rscala" % "2.5.0", // may not be necessary.
  "org.scalatest" %% "scalatest" % "3.0.0" % "test"
)

// Github sbt source dependencies.
//lazy val distribution = RootProject(uri("git://github.com/luiarthur/scala-distributions.git"))
//lazy val root = Project("root", file(".")).dependsOn(distribution)
